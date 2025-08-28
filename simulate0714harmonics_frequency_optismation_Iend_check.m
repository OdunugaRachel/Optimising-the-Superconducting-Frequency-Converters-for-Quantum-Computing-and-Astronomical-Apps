 
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
    
twpa.len = 110e-6*879*1/5;   % length of the the device, I believe this is fixed
% where does this length come from and whta does each factor represent?
twpa.betanl = 1;   %nonlinearity parameter, affects the strength of the nonlinear effects 
                   % %Josesphon Junction = 0.5

%% Define Simulation Sweep
% Frequencies over which to calculate gain. Defined here to set harmonic limit.
fcalc = 0.1e9:0.01e9:12.5e9;  % pump frequencies

%% Pick Modes


% Determine the highest possible harmonic that can be simulated for the *entire*
% pump frequency sweep, based on the available frequency range in the loaded data.
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

% array to store the complex current amplitude at the device output for each frequency
Iend = zeros(length(fcalc), length(twpa.modes)); 

for ii = 1:length(fcalc)
    twpa.pumpF = fcalc(ii);
%    If using dual pump instead use:
%    wn = twpa.modes(:,1)*twpa.pumpF + twpa.modes(:,2)*fcalc(ii) + twpa.modes(:,3)*twpa.pumpF2;

    % solveCME simulates the evolution of the current, including all linear
    % (loss, dispersion) and nonlinear effects. The output Y is the
    % physically correct complex current amplitude at each point z.
    Y = solveCME(fcalc(ii),zcalc,twpa);

    % Store the complex current at the end of the device (z = twpa.len)
    Iend(ii, :) = Y(end, :);

    disp(ii/length(fcalc))
end



%%  ANALYSIS & PLOTTING


position_plot_freq_GHz = 10.5; % Set the desired frequency in GHz for the position plot.

% --- Correct Power Calculation ---
% The output power relative to the input pump is calculated from the
% output current (Iend) and the input pump current (twpa.Ip).
output_power_dB = 20*log10(abs(Iend ./ twpa.Ip));

% --- Plot 1: Output Power vs. Frequency ---
disp('Generating Output Power vs. Frequency plot...');
f1 = figure('Name', 'Harmonic Output Power vs. Frequency');
hold on; grid on; box on;

plotLegend = {};
colors = lines(ceil(maxHarmonic/2)); % Use a standard color map

for i = 1:ceil(maxHarmonic/2)
    harmonic_num = 1 + 2*(i-1);
    plot(fcalc./1e9, smooth(output_power_dB(:,i), 1), 'LineWidth', 2, 'Color', colors(i,:));
    plotLegend{end+1} = sprintf('%dp', harmonic_num);
end

legend(plotLegend, 'Location', 'best');
xlim([fcalc(1) fcalc(end)]./1e9);
ylim([-40 10]);
xlabel('Pump Frequency (GHz)'); 
ylabel('Efficiency (dB)');
title('Harmonic Output Power vs. Pump Frequency');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Position', [100 100 1200 800]);
drawnow;

%% --- Plot 2: Power vs. Position at Specified Frequency ---

disp('Generating Power vs. Position plot...');

% Find the index of the frequency in fcalc that is closest to the desired one.
[~, freq_idx_for_plot] = min(abs(fcalc - position_plot_freq_GHz * 1e9));
actual_plot_freq_ghz = fcalc(freq_idx_for_plot) / 1e9;

% Re-run the simulation for the single specified frequency to get data vs. z
fprintf('Re-running simulation at %.3f GHz...\n', actual_plot_freq_ghz);
twpa.pumpF = fcalc(freq_idx_for_plot);
Y_position = solveCME(fcalc(freq_idx_for_plot), zcalc, twpa);

% Correctly calculate power vs position
power_vs_position_dB = 20*log10(abs(Y_position ./ twpa.Ip));

f2 = figure('Name', 'Power vs. Position');
hold on; grid on; box on;

for i = 1:ceil(maxHarmonic/2)
    plot(zcalc./twpa.len, smooth(power_vs_position_dB(:,i), 1), 'LineWidth', 2, 'Color', colors(i,:));
end

legend(plotLegend, 'Location', 'best');
xlabel('Normalized Position');
ylabel('Efficiency (dB)');
title(sprintf('Power vs. Position at %.2f GHz', actual_plot_freq_ghz));

set(gca, 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Position', [1350 100 1200 800]);
drawnow;


%% --- Save Plots ---
disp('--- Saving plots ---');
output_dir = 'Combined_Simulation_Results';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

saveas(f1, fullfile(output_dir, 'Output_Power_vs_Frequency.png'));
if exist('f2', 'var') && isvalid(f2)
    saveas(f2, fullfile(output_dir, 'Power_vs_Position.png'));
end
fprintf('Plots saved to the ''%s'' directory.\n', output_dir);

toc
