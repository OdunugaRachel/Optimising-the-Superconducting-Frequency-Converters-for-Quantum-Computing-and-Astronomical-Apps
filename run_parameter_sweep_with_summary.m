function run_parameter_sweep_with_summary()
% This script sweeps a geometry parameter and generates a summary plot.
%
%   This script systematically varies a single geometry parameter (e.g., 
%   modperiod, modamp, or unit_cells), runs the full simulation for each 
%   value to find the optimal pump frequency, and then creates a final 
%   summary plot visualizing the results and also generates plots for each simulation.
%
%   Instructions:
%   1. Choose ONE parameter to sweep by uncommenting the desired block below.
%   2. Set the `sweep_values` array to the values you want to test.
%   3. Set the `base_params` for the other two parameters.
%   4. Set the `pump_ratio` and `harmonic_to_optimize_for`.
%   5. Run the script.

close all;

%% Setting up the sweep

% --- Choose ONE parameter to sweep by uncommenting the desired block ---

% % 1. Sweep Modulation Period
% param_to_sweep = 'modperiod';
% sweep_values = linspace(10, 100, 40);
% base_params.modamp = 8.28;
% base_params.unit_cells = 300;

% % 2. Sweep Modulation Amplitude
% param_to_sweep = 'modamp';
% sweep_values = linspace(9, 13, 20);
% base_params.modperiod = 45;
% base_params.unit_cells = 813;

% 3. Sweep Unit Cells
param_to_sweep = 'unit_cells';
sweep_values = linspace(500, 900, 40);
base_params.modperiod = 120;
base_params.modamp = 3.9;

% Set common parameters for the sweep
pump_ratio = 0.13; % A single, constant pump power ratio for this sweep
harmonic_to_optimize_for = 3; % Which harmonic to optimize for (e.g., 3, 5)

%%  Sweeping parameter and running simulations

output_folder_name = ['Sweep_vs_', param_to_sweep];

% Pre-allocate arrays to store results for the final summary plot
optimal_frequencies = zeros(size(sweep_values));
max_powers = zeros(size(sweep_values));

for i = 1:length(sweep_values)
    geometry_params = base_params;
    geometry_params.(param_to_sweep) = sweep_values(i);
    fprintf('\n\n===== Running simulation for %s = %.2f =====\n', param_to_sweep, sweep_values(i));
    
    [opt_freq, max_pwr] = run_full_twpa_simulation_optimal_frequency(geometry_params, pump_ratio, output_folder_name, harmonic_to_optimize_for);
    
    % Store the results
    optimal_frequencies(i) = opt_freq;
    max_powers(i) = max_pwr;
end

%% SUMMARY PLOT

valid_indices = ~isnan(optimal_frequencies);
if ~any(valid_indices), warning('No valid optimal frequencies found. Skipping summary plot.'); return; end

f_summary = figure;
colororder({'b','r'});
scatter(sweep_values(valid_indices), optimal_frequencies(valid_indices), 75, 'filled', 'DisplayName', 'Optimal Frequency');
hold on; 
plot(sweep_values(valid_indices), optimal_frequencies(valid_indices), '--b', 'HandleVisibility', 'off'); 
grid on;
xlabel(['Swept Parameter: ', strrep(param_to_sweep, '_', ' ')]);
ylabel(sprintf('Optimal Pump Frequency for %dp (GHz)', harmonic_to_optimize_for));
title(['Optimization Summary vs. ', strrep(param_to_sweep, '_', ' ')]);
yyaxis right;
plot(sweep_values(valid_indices), max_powers(valid_indices), '-s', 'LineWidth', 2, 'MarkerFaceColor', 'auto', 'DisplayName', 'Max Power');
ylabel(sprintf('Max Power for %dp (dB)', harmonic_to_optimize_for));
legend('Location', 'best');
set(gca,'FontSize',12,'FontWeight','bold');

summary_filename = sprintf('Summary_vs_%s_for_%dp.png', param_to_sweep, harmonic_to_optimize_for);

% change directory to where you want to save the summary plot
saveas(f_summary, fullfile('Broader Frequency Range Successes', output_folder_name, summary_filename));
close(f_summary);

disp('Parameter sweep finished. Summary plot saved.');
end

