% This script simulates the harmonic generation in a TWPA (Traveling Wave Parametric Amplifier)
% for harmonic that can be supported by the frequency range of the loaded data.
% Optimizes the pump frequency for maximum output at a chosen harmonic.

clear

%% Initialize TWPA
tic
% load paramp_NbTiN_JPL_mK_short879_259_281_x0.875_taper_highf.mat
load 0714TWPaX.mat
%  containing previously calculated frequency (f), propagation constant (kperm), and transmission (S21)

% data = importdata("C:\Users\klimovich\Downloads\E08_06_750Cells_CAngDeg.csv");

twpa = createTWPA;

% twpa.fsim = data(:,1).*1e9;
% twpa.ksim = -data(:,2)/360*2*pi/twpaLEN;

% data = importdata("C:\Users\klimovich\Downloads\E08_06_750Cells_S21.csv");

% twpa.gsim = -data(:,3)/twpaLEN;
% twpa.gsim = twpa.gsim - twpa.gsim;


twpa.fsim = f; % array of frequencies
twpa.ksim = kperm;  % propagation constants
twpa.gsim = -log(abs(transpose(S21)));  % gain/loss per unit length (dB/m)

% Sanitize NaNs
% cleaning up NaNs in the frequency, k, and g arrays
twpa.ksim(isnan(twpa.ksim)) = max(twpa.ksim); % replacing nan values with max value
twpa.gsim(isnan(twpa.gsim)) = -100;

% twpa.pumpF = 6.3e9;
% twpa.pumpF2 = 0e9;
twpa.Istar = 4.5*1000;  % scaling factor for the current
%twpa.Ip = twpa.Istar*0.08; % pump current
twpa.Ip = twpa.Istar*0.13;
% processing time increases a lot as we set the pump current higher
twpa.Idc = twpa.Istar*0.0; % dc Bias current which is set to zero

len_scale = 1; % scaling factor for the length
    
twpa.len = (110e-6*879*1)/len_scale; % length of the the device, I believe this is fixed
% where does this length come from and whta does each factor represent?
twpa.betanl = 1;   %nonlinearity parameter, affects the strength of the nonlinear effects 
                   % %Josesphon Junction = 0.5

%% Define Simulation Sweep
% Frequencies over which to calculate gain. Defined here to set harmonic limit.
fcalc = 0.1e9:0.01e9:5.1e9;  % pump frequencies

%% Pick Modes


% Determine the highest possible harmonic that can be simulated for the *entire*
% pump frequency array, based on the available frequency range in the loaded data.
potential_harmonics = 9:-2:1; % Check from highest to lowest
max_pump_freq = max(fcalc);
max_available_freq = max(twpa.fsim);
maxHarmonic = 1; % Default to 1 (fundamental only)
 
for h = potential_harmonics
    if h * max_pump_freq <= max_available_freq
        maxHarmonic = h; % Set the new limit
        fprintf('Harmonic limit automatically set to %d based on f_calc range and loaded data.\n', maxHarmonic);
        break; % Found the highest possible harmonic, so we can stop.
    end
end

twpa.modes = [1 0]; % set up for harmonic modes (e.g. 1p, 3p, ..., 9p) - all odd harmonics up till the max
for i=3:2:maxHarmonic
    twpa.modes = cat(1, twpa.modes, [i 0]);
end

disp(twpa.modes)

% Alternative for dual-pump: with each mode being a row of 
% [(# of pump1 photons) (# of signal photons) (# of pump 2 photons)]
% so the example below is "p1, s, p2, 2p1 - s, 2p2 - s
% twpa.modes = [1 0 0; 0 1 0; 0 0 1; 2 -1 0; 0 -1 2];
%
twpa.I0 = zeros(length(twpa.modes),1);
twpa.I0(1) = twpa.Ip;

% twpa.I0(2) = twpa.Ip.*1e-3;


%% Compute Gain

% Positions to sample
zcalc = 0:0.0001:twpa.len; % positions along the TWPA length

g = zeros(length(fcalc), length(zcalc), length(twpa.modes));  % gain array, creating empty array to store the gain values
Iend = zeros(length(fcalc), length(twpa.modes)); % To store the complex current at the end of the device

for ii = 1:length(fcalc)
    wn = twpa.modes(:,1)*fcalc(ii);
    twpa.pumpF = fcalc(ii);
    S21_linear_prop = exp((-twpa.g(wn.') + 1i.*twpa.k(wn.')).*twpa.len);

    Y = solveCME(fcalc(ii),zcalc,twpa);

    g(ii,:,:) = 20*log10(abs(Y(:,:).*S21_linear_prop./twpa.I0(1)));  
    disp(ii/length(fcalc))
end

%% Plot Results
% close all

harmonic_to_optimize_for = 3; % which harmonic to optimize the pump frequency for

f3 = figure(3);
hold all
plotLegend = {};
for i=1:ceil(maxHarmonic/2)
    plot(fcalc./1e9,smooth(g(:,end,i),1),'Linewidth',2)
    plotLegend = [plotLegend, {[num2str(1 + 2*(i-1)),'p']}];
end
legend(plotLegend, 'Location', 'best')
grid on
xlim([fcalc(1) fcalc(end)]./1e9)
ylim([-30 5])
xlabel('Frequency (GHz) ')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
title('Output Power vs Frequency')
drawnow


% Finding the pump frequency that produces the maximum power for the desired harmonic.

% Calculate the index corresponding to the desired harmonic in the 'modes' dimension
mode_to_optimize_idx = (harmonic_to_optimize_for + 1) / 2;

% Check if this harmonic was actually simulated or is an invalid choice (e.g., even number)
if mode_to_optimize_idx > size(g, 3) || mod(harmonic_to_optimize_for, 2) == 0
    warning('Harmonic %d was not simulated or is invalid. Plotting for the last frequency instead.', harmonic_to_optimize_for);
    optimal_freq_idx = length(fcalc);
    optimal_pump_freq_GHz = fcalc(optimal_freq_idx) / 1e9;
    fprintf('Plotting for last simulated pump frequency: %.3f GHz\n', optimal_pump_freq_GHz);
else
    % Find the row (frequency index) that has the maximum power for the chosen harmonic at the output
    output_power_of_harmonic = g(:, end, mode_to_optimize_idx);
    [max_power_dB, optimal_freq_idx] = max(output_power_of_harmonic);
    optimal_pump_freq_GHz = fcalc(optimal_freq_idx) / 1e9;
    
    fprintf('Optimized for %d-th harmonic.\n', harmonic_to_optimize_for);
    fprintf('Maximum %dp power: %.2f dB\n', harmonic_to_optimize_for, max_power_dB);
    fprintf('This occurs at a pump frequency of: %.3f GHz\n', optimal_pump_freq_GHz);
end

f4 = figure(4);
cmap = colormap('turbo');
hold all
for i=1:ceil(maxHarmonic/2)
    plot(zcalc./twpa.len,smooth(g(optimal_freq_idx,:,i),1),'Linewidth',2)
end
legend(plotLegend, 'Location', 'best')
grid on
xlabel('Position')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
title(sprintf('Output Power vs Position at Optimal Pump Freq for %dp (%.2f GHz)', harmonic_to_optimize_for, optimal_pump_freq_GHz));
drawnow

%% 

saveas(f3, 'Frequency_Power.png')
saveas(f4, 'BestFrequency_Harmonics.png')

toc

