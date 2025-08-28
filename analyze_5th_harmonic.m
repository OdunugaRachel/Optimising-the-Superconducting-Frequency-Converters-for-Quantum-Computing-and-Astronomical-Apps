% This script analyzes the 5th harmonic generation in a TWPA (Traveling Wave Parametric Amplifier).
% It includes the calculation of phase mismatch (delta beta) and efficiency of the 5th harmonic generation
% and plots delta beta against efficiency, marking the point of maximum efficiency.

clear;
close all;
tic;

% I have an inkling that a bunch of this code is duplicated from the other files and a bit redundant but it still works so I'm leaving it for now
%% 1. INITIALIZE TWPA AND SIMULATION PARAMETERS

% Load pre-calculated network properties
load 0714TWPaX.mat;

% Create and configure the TWPA object
twpa = createTWPA;
twpa.fsim = f;
twpa.ksim = kperm;
twpa.gsim = -log(abs(transpose(S21)));

% Sanitize NaNs from the input data
twpa.ksim(isnan(twpa.ksim)) = max(twpa.ksim);
twpa.gsim(isnan(twpa.gsim)) = -100;

% Set TWPA physical and pump parameters
twpa.Istar = 4.5*1000;  % Scaling factor for the current
twpa.Ip = twpa.Istar * 0.13; % Pump current (adjust as needed)
twpa.Idc = twpa.Istar * 0.0; % DC Bias current
twpa.len = len_meters; % Use length from the loaded .mat file for consistency
twpa.betanl = 1;       % Nonlinearity parameter

%% 2. DEFINE MODES AND SIMULATION SWEEP

% Define which harmonic to analyze. This makes the script easily adaptable.
harmonic_to_analyze = 5;

% Automatically determine the maximum harmonic possible based on data range
fcalc = 0.1e9:0.01e9:5.1e9;  % Pump frequencies to sweep
max_pump_freq = max(fcalc);
max_available_freq = max(twpa.fsim);
potential_harmonics = 9:-2:1;
maxHarmonic = 1;
for h = potential_harmonics
    if h * max_pump_freq <= max_available_freq
        maxHarmonic = h;
        fprintf('Harmonic limit automatically set to %d based on f_calc range and loaded data.\n', maxHarmonic);
        break;
    end
end

if harmonic_to_analyze > maxHarmonic
    error('The harmonic to analyze (%d) is higher than the max possible harmonic (%d) for the given frequency sweep.', harmonic_to_analyze, maxHarmonic);
end

% Set up modes for simulation (1p, 3p, 5p, etc., up to maxHarmonic)
twpa.modes = [1 0];
for i = 3:2:maxHarmonic
    twpa.modes = cat(1, twpa.modes, [i 0]);
end
disp('Simulating modes:');
disp(twpa.modes);

% Set initial current conditions (only pump is on at z=0)
twpa.I0 = zeros(length(twpa.modes), 1);
twpa.I0(1) = twpa.Ip;

%%  3. COMPUTE HARMONIC POWER (MAIN SIMULATION LOOP)

% Positions to sample along the device
zcalc = 0:0.0001:twpa.len;

% Pre-allocate gain array
g = zeros(length(fcalc), length(zcalc), length(twpa.modes));

%  starting a parallel pool to speed up the simulation
if isempty(gcp('nocreate'))
    disp('Starting parallel pool...');
    parpool;
end

parfor ii = 1:length(fcalc)
    twpa_local = twpa; % Create a local copy for parfor
    twpa_local.pumpF = fcalc(ii);
    wn = twpa_local.modes(:,1) * fcalc(ii);
    
    % Calculate S21 for the current set of frequencies
    S21_val = exp((-twpa_local.g(wn.') + 1i.*twpa_local.k(wn.')).*twpa_local.len);

    Y = solveCME(fcalc(ii), zcalc, twpa_local);

    % Calculate power in dB relative to the input pump power
    g(ii,:,:) = 20*log10(abs(Y(:,:).*S21_val./twpa_local.I0(1)));
    
    if mod(ii, 20) == 0
        fprintf('Progress: %.1f%%\n', (ii/length(fcalc))*100);
    end
end
disp('--- Main simulation finished ---');

%%  4. ANALYSIS PART A: FIND OPTIMAL PUMP FREQUENCY

% Calculate the index for the harmonic we want to analyze
% e.g., for 5p, modes are [1 0; 3 0; 5 0; ...], so 5p is the 3rd mode.
mode_idx = (harmonic_to_analyze + 1) / 2;

% Find the pump frequency that maximizes the output power of the target harmonic
output_power_of_harmonic = g(:, end, mode_idx);
[max_power_dB, optimal_freq_idx] = max(output_power_of_harmonic);
optimal_pump_freq_GHz = fcalc(optimal_freq_idx) / 1e9;

% Display results
fprintf('\n--- Optimal Pump Frequency for %dp Harmonic ---\n', harmonic_to_analyze);
fprintf('  - Maximum %dp Power: %.2f dB\n', harmonic_to_analyze, max_power_dB);
fprintf('  - Optimal Pump Frequency: %.3f GHz\n', optimal_pump_freq_GHz);

% --- Plotting: Power vs. Position at Optimal Frequency ---
f1 = figure(1);
hold all;
plotLegend = {};
for i = 1:ceil(maxHarmonic/2)
    harmonic_num = 1 + 2*(i-1);
    plot(zcalc./twpa.len, smooth(g(optimal_freq_idx,:,i),1), 'Linewidth', 2);
    plotLegend = [plotLegend, {[num2str(harmonic_num),'p']}];
end
legend(plotLegend);
grid on;
xlabel('Normalized Position (z/L)');
ylabel('Power (dB)');
title(sprintf('Harmonic Power vs. Position at Optimal Pump Freq for %dp (%.2f GHz)', harmonic_to_analyze, optimal_pump_freq_GHz));
set(gca,'FontSize',14,'FontWeight','bold');
set(gcf,'Position',[100 100 1200 800]);
drawnow;

%%  5. ANALYSIS PART B: CONVERSION EFFICIENCY VS. PHASE MISMATCH

% Set up arrays for delta beta and efficiency
delta_beta = zeros(size(fcalc));
efficiency = zeros(size(fcalc));

for ii = 1:length(fcalc)
    % Interpolate to find k at the pump and harmonic frequencies
    k_pump  = interp1(twpa.fsim, twpa.ksim, fcalc(ii));
    k_Nharm = interp1(twpa.fsim, twpa.ksim, harmonic_to_analyze * fcalc(ii));
    
    % Calculate phase mismatch: Δβ = k(Nω) - N*k(ω)
    delta_beta(ii) = k_Nharm - harmonic_to_analyze * k_pump;
    
    % Calculate efficiency: P_out(Nω) / P_in(ω)
    % g is power in dB relative to P_in(ω), so P_out(Nω)/P_in(ω) = 10^(g/10)
    efficiency(ii) = 10^(g(ii, end, mode_idx) / 10);
end

% Find the maximum efficiency point
[max_efficiency, max_idx] = max(efficiency);
max_delta_beta = delta_beta(max_idx);

% Display results
fprintf('\n--- Maximum Efficiency Point for %dp Harmonic ---\n', harmonic_to_analyze);
fprintf('  - Max Conversion Efficiency: %.4f (linear)\n', max_efficiency);
fprintf('  - Phase Mismatch (Δβ) at Max Efficiency: %.2f rad/m\n', max_delta_beta);

% --- Plotting: Efficiency vs. Phase Mismatch ---
f2 = figure(2);
scatter(delta_beta, efficiency, 'filled', 'DisplayName', 'Data');
hold on;
plot(max_delta_beta, max_efficiency, 'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Max Efficiency');
text_label = sprintf('Max Eff: %.3f\n at \\Delta\\beta = %.2f', max_efficiency, max_delta_beta);
text(max_delta_beta, max_efficiency, text_label, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12, 'Color', 'r');
grid on;
xlabel('Phase Mismatch, \Delta\beta (rad/m)');
ylabel(sprintf('%dth Harmonic Conversion Efficiency (Linear Power)', harmonic_to_analyze));
title(sprintf('\\Delta\\beta vs. %dth Harmonic Efficiency', harmonic_to_analyze));
legend('show', 'Location', 'best');
set(gca,'FontSize',14,'FontWeight','bold');
set(gcf,'Position',[1350 100 1200 800]);
drawnow;

%%   6. SAVE RESULTS

% Create a directory for the output if it doesn't exist
output_dir = '5th_Harmonic_Analysis_Results';
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end

% Save the figures
saveas(f1, fullfile(output_dir, 'OptimalFrequency_Power_vs_Position_5p.png'));
saveas(f2, fullfile(output_dir, 'Efficiency_vs_DeltaBeta_5p.png'));

fprintf('Plots saved to the ''%s'' directory.\n', output_dir);

toc;

