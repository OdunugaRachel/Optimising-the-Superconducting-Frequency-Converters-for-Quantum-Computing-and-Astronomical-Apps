function run_full_twpa_simulation(geometry_params, pump_ratio, output_folder_name)
%   This file combines network calculation and harmonic simulation into a single file.
%
%   This script first calculates the physical properties (dispersion, S21) of a
%   TWPA based on its geometry, and then immediately uses those properties to
%   simulate harmonic generation for a range of pump powers.
%
%   INPUTS:
%   geometry_params: A struct with fields 'modperiod', 'modamp', 'unit_cells'.
%   pump_ratio: A scalar value for the pump power ratio (Ip/Istar).
%   output_folder_name: The name of the subfolder to save plots in.
%
%   OUTPUTS:
%   Plot files are saved in the specified output folder.

% clear; 
% close all;
tic;

%%   PART 1: NETWORK CALCULATION (from network_calc_adaptive_freq.m)

% Define geometry and modulation parameters.
g.finger.l = 20;
g.finger.w = 0.25;
g.finger.p = 0.87;
g.ms.w = g.finger.w;
g.ms.h = 0.15;
g.ms.eps = 10.3;
g.ms.eps_upper = 10.4;
g.ms.t1 = 0.030;
g.ms.t2 = 0.200;

% This comes from the sweepin script that gives the geometry_params struct
g.finger.modperiod = geometry_params.modperiod;
g.finger.modamp = geometry_params.modamp;
unit_cells = geometry_params.unit_cells;
L = g.finger.modperiod * unit_cells; % Total length in um

% Define physical constants and simulation parameters
f_initial = linspace(.01, 50.01, 60000) * 1e9;
rhon = 480.*1e-8;
sn = 1 / rhon;
ev = 1.6e-19;
mu0 = pi*4e-7;
kb = 1.38e-23;
T = 0.01;
kbt = kb * T / ev;
Tc = 14.5;
gap = 0.5 * 3.5 * Tc * kb / ev;
ms = g.ms;

% Run adaptive frequency grid calculation
phase_jump_threshold = 0.85 * pi; % Threshold for phase jumps, change if needed
[f, S21, ~, ~, ~, ~] = refine_freq_grid(f_initial, g, ms, gap, kbt, mu0, sn, L, phase_jump_threshold);

uth = -1*unwrap(angle(S21(:)));
for fj = 2:length(f)
    while uth(fj)<uth(fj-1)
        uth(fj:end) = uth(fj:end) + 2*pi;
    end
end
slope = (uth(2)-uth(1))./(f(2)-f(1));
intercept = uth(1) - slope*f(1);
uth = uth - intercept;
len_meters = L / 1e6; % Convert length to meters
kperm = uth / len_meters;

%%  PART 2: HARMONIC SIMULATION (from simulate0714harmonics_series.m)

% Initialize TWPA struct
twpa = createTWPA;

% Populate with data from Part 1 (instead of loading from file)
twpa.fsim = f;
twpa.ksim = kperm;
twpa.gsim = -log(abs(transpose(S21)));
twpa.len = len_meters; % Use the calculated length for consistency

% Sanitize NaNs
twpa.ksim(isnan(twpa.ksim)) = max(twpa.ksim);
twpa.gsim(isnan(twpa.gsim)) = -100;

% Define simulation parameters
twpa.Istar = 4.5*1000;
twpa.betanl = 1;

% Pick Modes
maxHarmonic = 9;
twpa.modes = [1 0];
for i=3:2:maxHarmonic
    twpa.modes = cat(1, twpa.modes, [i 0]);
end

disp(twpa.modes);
twpa.I0 = zeros(length(twpa.modes),1);

% starting a parallel pool so we can speed up the simulations
if isempty(gcp('nocreate'))
    disp('Starting parallel pool...');
    parpool;
end
fcalc = 0.1e9:0.01e9:5.1e9;
zcalc = 0:0.0001:twpa.len;
results = struct(); % Initialize results as a structure


current_ratio = pump_ratio;
twpa_iter = twpa;
twpa_iter.Ip = current_ratio * twpa_iter.Istar;
twpa_iter.I0(1) = twpa_iter.Ip;
g_data = zeros(length(fcalc), length(zcalc), length(twpa.modes));
parfor ii = 1:length(fcalc)
    twpa_local = twpa_iter;
    twpa_local.pumpF = fcalc(ii);
    wn = twpa_local.modes(:,1) * fcalc(ii);
    S21_val = exp((-twpa_local.g(wn.') + 1i.*twpa_local.k(wn.')).*twpa_local.len);
    Y = solveCME(fcalc(ii), zcalc, twpa_local);
    g_data(ii,:,:) = 20*log10(abs(Y(:,:).*S21_val./twpa_local.I0(1)));

end

% Store results directly, without the cell array wrapper
results.g = g_data;
results.ratio = current_ratio;

% plotting and saving results

% change directory to where you want to save the plots
output_dir = fullfile('Broader Frequency Range Successes', output_folder_name);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Get data from our stored results
g_data = results.g;
current_ratio = results.ratio;

% create figure, hidden to avoid pop-up during batch runs
f_handle = figure('Visible', 'off', 'Position', [100 100 1200 1000]);

% Subplot 1 (Top-Left): S21 vs Frequency
subplot(2, 2, 1);
plot(f./1e9, abs(S21), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
grid on;
title('Transmission (S_{21})');
xlabel('Frequency (GHz)');
ylabel('|S_{21}|');
set(gca,'FontSize',10,'FontWeight','bold');

% Subplot 2: Dispersion vs Frequency
subplot(2, 2, 2);
plot(f./1e9, kperm, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
grid on;
title('Dispersion');
xlabel('Frequency (GHz)');
ylabel('k_{perm} (rad/m)');
set(gca,'FontSize',10,'FontWeight','bold');

% Subplot 3: Output Power vs Frequency
subplot(2, 2, 3);
hold all;
plotLegend = {};
for h_idx = 1:ceil(maxHarmonic/2)
    plot(fcalc./1e9, smooth(g_data(:,end,h_idx),1), 'Linewidth', 2);
    plotLegend = [plotLegend, {[num2str(1 + 2*(h_idx-1)),'p']}];
end
    legend(plotLegend, 'Location', 'best');
    grid on; xlim([fcalc(1) fcalc(end)]./1e9); ylim([-30 5]); xlabel('Frequency (GHz)'); ylabel('Output Power (dB)'); title('Harmonic Output Power vs. Frequency'); set(gca,'FontSize',10,'FontWeight','bold');

% Subplot 4: Output Power vs Position
subplot(2, 2, 4);
hold all;
index = length(fcalc);
for h_idx = 1:ceil(maxHarmonic/2)
    plot(zcalc./twpa.len, smooth(g_data(index,:,h_idx),1), 'Linewidth', 2);
end
legend(plotLegend, 'Location', 'best'); grid on; xlabel('Position (Normalized)'); ylabel('Output Power (dB)'); title('Harmonic Output Power vs. Position'); set(gca,'FontSize',10,'FontWeight','bold');

main_title = sprintf('Full Sim: modperiod=%.1f, modamp=%.1f, cells=%d, ratio=%.3f', g.finger.modperiod, g.finger.modamp, unit_cells, current_ratio);
sgtitle(f_handle, main_title);
output_filename = sprintf('Full_Sim_modperiod_%.1f_modamp_%.1f_cells_%d_ratio_%.3f.png', g.finger.modperiod, g.finger.modamp, unit_cells, current_ratio);
saveas(f_handle, fullfile(output_dir, output_filename));
close(f_handle);

disp('All plots saved successfully.');
toc;

end
