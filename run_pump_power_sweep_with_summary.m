function run_pump_power_sweep_with_summary()

close all;

%% 

param_to_sweep = 'pump_ratio';
sweep_values = linspace(0.05, 0.15, 20); % Example: sweep from 5% to 15% of Istar

geometry_params.modperiod = 90;
geometry_params.modamp = 7;
geometry_params.unit_cells = 830;


harmonic_to_optimize_for = 3; % Which harmonic to optimize for (e.g., 3, 5)

%% 
output_folder_name = ['Sweep_vs_', param_to_sweep];

optimal_frequencies = zeros(size(sweep_values));
max_powers = zeros(size(sweep_values));

for i = 1:length(sweep_values)
    current_pump_ratio = sweep_values(i);
    fprintf('\n\n===== Running simulation for %s = %.3f =====\n', param_to_sweep, current_pump_ratio);
    
    [opt_freq, max_pwr] = run_full_twpa_simulation_optimal_frequency(geometry_params, current_pump_ratio, output_folder_name, harmonic_to_optimize_for);
    
    % Store the results
    optimal_frequencies(i) = opt_freq;
    max_powers(i) = max_pwr;
end

%% Summary Plot
disp('--- Generating Summary Plot ---');

valid_indices = ~isnan(optimal_frequencies);
if ~any(valid_indices), warning('No valid optimal frequencies found. Skipping summary plot.'); return; end

f_summary = figure;
colororder({'b','r'});
scatter(sweep_values(valid_indices), optimal_frequencies(valid_indices), 75, 'filled', 'DisplayName', 'Optimal Frequency');
hold on; plot(sweep_values(valid_indices), optimal_frequencies(valid_indices), '--b', 'HandleVisibility', 'off'); grid on;
xlabel('Pump Power Ratio (Ip/Istar)');
ylabel(sprintf('Optimal Pump Frequency for %dp (GHz)', harmonic_to_optimize_for));
title(['Optimization Summary vs. Pump Power Ratio']);
yyaxis right;
plot(sweep_values(valid_indices), max_powers(valid_indices), '-s', 'LineWidth', 2, 'MarkerFaceColor', 'auto', 'DisplayName', 'Max Power');
ylabel(sprintf('Max Power for %dp (dB)', harmonic_to_optimize_for));
legend('Location', 'best');
set(gca,'FontSize',12,'FontWeight','bold');

summary_filename = sprintf('Summary_vs_%s_for_%dp.png', param_to_sweep, harmonic_to_optimize_for);
saveas(f_summary, fullfile('Secure Third Harmonic Generation', output_folder_name, summary_filename));
close(f_summary);

disp('Parameter sweep finished. Summary plot saved.');
end

