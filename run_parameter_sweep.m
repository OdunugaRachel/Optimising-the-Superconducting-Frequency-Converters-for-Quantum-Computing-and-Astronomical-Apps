function run_parameter_sweep()
% RUN_PARAMETER_SWEEP Sweeps a single geometry parameter and runs the full simulation for each value.
%
%   This script allows you to easily test the effect of changing one
%   parameter (e.g., modperiod, modamp, or unit_cells) while keeping others
%   constant.
%
%   Instructions:
%   1. Choose ONE parameter to sweep by uncommenting the desired block below.
%   2. Set the `sweep_values` array to the values you want to test.
%   3. Set the `base_params` for the other two parameters.
%   4. Set the `pump_ratio` to a constant value for the entire sweep.
%   5. Run the script.

close all;

%% ========================================================================
%  SWEEP CONFIGURATION
%  ========================================================================

% --- Choose ONE parameter to sweep by uncommenting the desired block ---

% % 1. Sweep Modulation Period
% param_to_sweep = 'modperiod';
% sweep_values = linspace(70,100,10);
% base_params.modamp = 7;
% base_params.unit_cells = 830;

% 2. Sweep Modulation Amplitude
param_to_sweep = 'modamp';
sweep_values = linspace(4,8,10);
base_params.modperiod = 90;
base_params.unit_cells = 830;

% % 3. Sweep Unit Cells
% param_to_sweep = 'unit_cells';
% sweep_values = linspace(400, 1000 ,20);
% base_params.modperiod = 90;
% base_params.modamp = 7;

% --- Set common parameters for the sweep ---
pump_ratio = 0.13; % A single, constant pump power ratio for this sweep
%harmonic_to_optimize_for = 3; % Which harmonic to optimize for (e.g., 3, 5), when using optimal frequency simulation

%% ========================================================================
%  RUN SWEEP (No changes needed below this line)
%  ========================================================================

% Determine output folder name based on the parameter being swept
switch param_to_sweep
    case 'modperiod'
        output_folder_name = 'Different Modulation Period';
    case 'modamp'
        output_folder_name = 'Different Modulation Amplitude';
    case 'unit_cells'
        output_folder_name = 'Different Unit Cells';
    otherwise
        % A sensible default if a new parameter is added without a case
        output_folder_name = [param_to_sweep, '_Sweep'];
end

for i = 1:length(sweep_values)
    geometry_params = base_params;
    geometry_params.(param_to_sweep) = sweep_values(i);
    fprintf('\n\n===== Running simulation for %s = %.2f =====\n', param_to_sweep, sweep_values(i));
    run_full_twpa_simulation(geometry_params, pump_ratio, output_folder_name); % for set frequency
    % Uncomment the line below to run the simulation with optimal frequency
    %[~, ~] = run_full_twpa_simulation_optimal_frequency(geometry_params, pump_ratio, output_folder_name, harmonic_to_optimize_for)
disp('Parameter sweep finished.');
end
