clear 
%close all
tic;

load 0714TWPaX.mat

f_pump = linspace(0.01e9, 15e9, 30000);

target_delta_beta = -22.79;
%target_delta_beta_5th = -21.88;



%% ========================================================================
harmonics_to_check = [3, 5];  % Now check both 3rd and 5th harmonics

delta_beta = zeros(length(harmonics_to_check), length(f_pump));

disp('--- Calculating Delta Beta for each harmonic ---');

prev_data = 0;

k_at_pump_freq = interp1(f, kperm, f_pump, 'linear', 'extrap');

for i = 1:length(harmonics_to_check)
    N = harmonics_to_check(i);
    
    fprintf('Calculating for %dp harmonic...\n', N);

    k_at_N_harmonic_freq = interp1(f, kperm, N * f_pump, 'linear', 'extrap');

    if N == 3
        delta_beta(i, :) = k_at_N_harmonic_freq - 3 * k_at_pump_freq;
    else
        
        sum_intermediate_k = zeros(size(f_pump));
        intermediate_harmonics = 3:2:(N-2);
        
        for m = intermediate_harmonics
            k_intermediate = interp1(f, kperm, m * f_pump, 'linear', 'extrap');
            sum_intermediate_k = sum_intermediate_k + k_intermediate;
        end
        
        delta_beta(i, :) = k_at_N_harmonic_freq - 2 * k_at_pump_freq - sum_intermediate_k;
    end
end
%% Plotting the results

disp('--- Generating Plots ---');

f_plot = figure(10); % Selects figure 10 or creates it if it doesn't exist.
clf(f_plot); % Clears the figure content to ensure a fresh plot.
set(f_plot, 'Name', 'Delta Beta Analysis', 'Position', [200 200 1200 800]);
hold on; grid on; box on;

plotLegend = {};
colors = lines(length(harmonics_to_check)); % Get distinct colors

for i = 1:length(harmonics_to_check)
    N = harmonics_to_check(i);

    plot(f_pump / 1e9, delta_beta(i, :), 'LineWidth', 2, 'Color', colors(i, :), 'LineStyle', '-');
    plotLegend{end+1} = sprintf('%dp Harmonic', N);
end

% Plot the target line and update its legend entry
yline(target_delta_beta, 'k--', 'DisplayName', sprintf('\\Delta\\beta = %.2f (Target)', target_delta_beta));
plotLegend{end+1} = sprintf('Target (%.2f)', target_delta_beta);

title('Phase Mismatch (\\Delta\\beta) vs. Pump Frequency for 3rd and 5th Harmonics');
xlabel('Pump Frequency (GHz)');
ylabel('Phase Mismatch, \\Delta\\beta (rad/m)');
legend(plotLegend, 'Location', 'best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([min(f_pump) max(f_pump)]/1e9);

% Save the plot
output_dir = 'Delta_Beta_Analysis';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_filename = 'DeltaBeta_vs_Frequency_3rd_5th.png';
saveas(f_plot, fullfile(output_dir, output_filename));

fprintf('Plot saved to %s\\%s\n', output_dir, output_filename);

% Optionally, plot 5th vs 3rd harmonic phase mismatch directly
figure(11); clf; hold on; grid on; box on;
scatter(delta_beta(1,:), delta_beta(2,:), 5, 'b', 'filled');
% Add target lines for both axes
xline(target_delta_beta, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Target 3rd (%.2f)', target_delta_beta));
yline(target_delta_beta, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Target 5th (%.2f)', target_delta_beta));
xlabel('3rd Harmonic Phase Mismatch (rad/m)');
ylabel('5th Harmonic Phase Mismatch (rad/m)');
title('5th Harmonic vs 3rd Harmonic Phase Mismatch');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best');

saveas(gcf, fullfile(output_dir, 'DeltaBeta_5th_vs_3rd.png'));

toc;