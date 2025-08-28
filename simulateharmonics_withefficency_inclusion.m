
clear

tic
%% Initialize TWPA

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
twpa.Ip = twpa.Istar*0.12;
% processing time increases a lot as we set the pump current higher
twpa.Idc = twpa.Istar*0.0; % dc Bias current which is set to zero
    
twpa.len = 110e-6*879*1;   % length of the the device, I believe this is fixed
% where does this length come from and whta does each factor represent?
twpa.betanl = 1;   %nonlinearity parameter, affects the strength of the nonlinear effects 
                   % %Josesphon Junction = 0.5


%% Pick Modes

maxHarmonic = 9;  %max harmonic order to simulate

twpa.modes = [1 0]; % set up for harmonic modes (e.g. 1p, 3p, ..., 9p) - all odd harmonics up till the max
for i=3:2:maxHarmonic
    twpa.modes = cat(1, twpa.modes, [i 0]);
end

disp(twpa.modes)

% Alternative for dual-pump: with each mode being a row of 
% [(# of pump1 photons) (# of signal photons) (# of pump 2 photons)]
% so the example below is "p1, s, p2, 2p1 - s, 2p2 - s
% twpa.modes = [1 0 0; 0 1 0; 0 0 1; 2 -1 0; 0 -1 2];

twpa.I0 = zeros(length(twpa.modes),1);
twpa.I0(1) = twpa.Ip;
% twpa.I0(2) = twpa.Ip.*1e-3;

%% Compute Gain

% Frequencies over which to calculate gain
fcalc = 0.1e9:0.005e9:5.1e9;  % pump frequencies

% Positions to sample
zcalc = 0:0.0001:twpa.len; % positions along the TWPA length

g = zeros(length(fcalc), length(zcalc), length(twpa.modes));  % gain array, creating empty array to store the gain values
Iend = zeros(length(fcalc), length(zcalc), length(twpa.modes)); % end current 

for ii = 1:length(fcalc)
    wn = twpa.modes(:,1)*fcalc(ii);
    twpa.pumpF = fcalc(ii);
%    If using dual pump instead use:
%    wn = twpa.modes(:,1)*twpa.pumpF + twpa.modes(:,2)*fcalc(ii) + twpa.modes(:,3)*twpa.pumpF2;
    S21 = exp((-twpa.g(wn.') + 1i.*twpa.k(wn.')).*twpa.len);

    Y = solveCME(fcalc(ii),zcalc,twpa);

    g(ii,:,:) = 20*log10(abs(Y(:,:).*S21./twpa.I0(1)));  
    disp(ii/length(fcalc))
end

%% Plot Results

% Plot frequency vs output power (existing)
f3 = figure(3);
hold all
plotLegend = {};
for i=1:ceil(maxHarmonic/2)
    plot(fcalc./1e9,smooth(g(:,end,i),1),'Linewidth',2)
    plotLegend = [plotLegend, {[num2str(1 + 2*(i-1)),'p']}];
end
legend(plotLegend)
grid on
xlim([fcalc(1) fcalc(end)]./1e9)
ylim([-30 5])
xlabel('Frequency (GHz) ')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
drawnow

% Plot position vs output power (existing)
index = length(fcalc);
f4 = figure(4);
cmap = colormap('turbo');
hold all
for i=1:ceil(maxHarmonic/2)
    plot(zcalc./twpa.len,smooth(g(index,:,i),1),'Linewidth',2)
end 
legend(plotLegend)
grid on
xlabel('Position')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
drawnow



%% Delta beta vs efficiency plot 


% Set up arrays for delta beta and efficiency
delta_beta = zeros(size(fcalc));
efficiency = zeros(size(fcalc));
k_pump = zeros(size(fcalc));
k_3harm = zeros(size(fcalc));

for ii = 1:length(fcalc)
    % Calculate k(3*omega_p) and k(omega_p) using interpolation
    k_pump(ii) = interp1(twpa.fsim, twpa.ksim, fcalc(ii));
    k_3harm(ii) = interp1(twpa.fsim, twpa.ksim, 3*fcalc(ii));
    delta_beta(ii) = k_3harm(ii) - 3*k_pump(ii);
    % Efficiency: use output power of 3rd harmonic at end of device
    % (Assume 3rd harmonic is i=2 in g, since twpa.modes = [1 0; 3 0; ...])
    efficiency(ii) = 10^(g(ii, end, 2)./10); % Power (dB) of 3rd harmonic at output
end

% Find the maximum efficiency and its corresponding delta beta
[max_efficiency, max_idx] = max(efficiency);
max_delta_beta = delta_beta(max_idx);
max_k_pump = k_pump(max_idx);
max_k_3harm = k_3harm(max_idx);


% Display the results in the command window
fprintf('\n--- Maximum Efficiency Point ---\n');
fprintf('  - Max Efficiency: %.4f (linear)\n', max_efficiency);
fprintf('  - Delta Beta at Max Efficiency: %.2f rad/m\n', max_delta_beta);
fprintf('  - k(3*omega_p) at Max Efficiency: %.4f rad/m\n', max_k_3harm);
fprintf('  - k(omega_p) at Max Efficiency: %.4f rad/m\n', max_k_pump);

% Plot delta beta vs efficiency
f5 = figure(5);
scatter(delta_beta, efficiency, 'filled', 'DisplayName', 'Data')
hold on;

% Mark the maximum point with a red star
plot(max_delta_beta, max_efficiency, 'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Max Efficiency');

% Add a text label to the max point
text_label = sprintf('Max Eff: %.3f\n at \\Delta\\beta = %.2f', max_efficiency, max_delta_beta);
text(max_delta_beta, max_efficiency, text_label, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12, 'Color', 'r');

grid on
xlabel('Phase Mismatch, \Delta\beta (rad/m)')
ylabel('3rd Harmonic Conversion Efficiency (Linear Power)')
title('\Delta\beta vs Efficiency')
legend('show', 'Location', 'best')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000]) % Use a valid position vector
drawnow
saveas(f5, 'DeltaBeta_vs_Efficiency.png')

% %% Sum of squared amplitudes vs position (energy conservation along device)
% sum_sq_pos = zeros(length(zcalc),1);
% for jj = 1:length(zcalc)
%     amps = squeeze(g(index,jj,:)); % dB values for all harmonics at this position
%     amps_lin = 10.^(amps/20);      % convert dB to linear amplitude
%     sum_sq_pos(jj) = sum(amps_lin.^2); % sum of squares
% end
% 
% f_sum_pos = figure;
% plot(zcalc./twpa.len, sum_sq_pos, 'k-', 'LineWidth', 2)
% grid on
% xlabel('Position (normalized)')
% ylabel('Sum of squared amplitudes')
% title('Sum of squared amplitudes vs position (should be constant if energy conserved)')
% set(gca,'FontSize',16)
% set(gca,'FontWeight','bold')
% set(gcf,'Position',[1500 100 1500 1000])
% drawnow
% saveas(f_sum_pos, 'SumSquaredAmplitudes_vs_Position.png')
% 
% saveas(f3, 'Frequency_Power.png')
% saveas(f4, 'Harmonics.png')

toc
