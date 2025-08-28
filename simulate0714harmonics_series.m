
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

%% --- Simulation Configuration ---
single_frequency_override = []; % Set to a value in Hz (e.g., 4.5e9) to override the frequency sweep.
ratio_range = linspace(0.03, 0.08, 8); % Range of Ip/Istar ratios to simulate.
harmonic_to_optimize_for = 3; % Which harmonic to optimize the 'Power vs Position' plot for.

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

twpa.I0 = zeros(length(twpa.modes),1);

% twpa.I0(2) = twpa.Ip.*1e-3;
%% 1. --- COMPUTATION STAGE ---
% We will calculate all results first and store them, separating computation from plotting.

% Start a parallel pool if one is not already running. Requires Parallel Computing Toolbox.
if isempty(gcp('nocreate'))
    disp('Starting parallel pool...');
    parpool;
end

% Frequencies over which to calculate gain.
% The simulation will always run over the full frequency range.
fcalc = 0.1e9:0.001e9:5.1e9;

% Positions to sample
zcalc = 0:0.0001:twpa.len; % positions along the TWPA length

% Pre-allocate a cell array to store the results for each ratio
results = cell(1, length(ratio_range));

for r_idx = 1:length(ratio_range)
    current_ratio = ratio_range(r_idx);
    fprintf('Calculating for ratio: %.3f\n', current_ratio);

    % Create a temporary twpa struct for this iteration to ensure parfor safety
    twpa_iter = twpa;
    twpa_iter.Ip = current_ratio * twpa_iter.Istar;
    twpa_iter.I0(1) = twpa_iter.Ip;

    % Preallocate gain array for this ratio
    g = zeros(length(fcalc), length(zcalc), length(twpa.modes));

    % Use parfor for the most expensive loop. This distributes frequency calculations across CPU cores.
    parfor ii = 1:length(fcalc)
        % Create a local copy of the twpa object for this specific frequency iteration
        twpa_local = twpa_iter;
        twpa_local.pumpF = fcalc(ii);

        wn = twpa_local.modes(:,1) * fcalc(ii);
        S21_val = exp((-twpa_local.g(wn.') + 1i.*twpa_local.k(wn.')).*twpa_local.len);

        Y = solveCME(fcalc(ii), zcalc, twpa_local);

        g(ii,:,:) = 20*log10(abs(Y(:,:).*S21_val./twpa_local.I0(1)));
    end

    % Store the results for this ratio
    results{r_idx}.g = g;
    results{r_idx}.ratio = current_ratio;

end

%% 2. --- PLOTTING STAGE ---
% Now that all calculations are done, loop through the results and create the plots.
disp('All calculations finished. Now generating plots...');
output_dir = fullfile('Ninth Harmonic Search', 'Different Pump Powers');
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end
for r_idx = 1:length(ratio_range)
    % Get data for this ratio from our stored results
    g_data = results{r_idx}.g;
    current_ratio = results{r_idx}.ratio;
    
    % Create figure, but keep it hidden while plotting for speed
    f_handle = figure('Visible', 'off');
    
    % --- Subplot 1: Output Power vs Frequency (Always plotted) ---
    subplot(2, 1, 1);
    hold all;
    plotLegend = {};
    for h_idx = 1:ceil(maxHarmonic/2)
        plot(fcalc./1e9, smooth(g_data(:,end,h_idx),1), 'Linewidth', 2);
        plotLegend{end+1} = [num2str(1 + 2*(h_idx-1)),'p'];
    end
    legend(plotLegend);
    grid on;
    xlim([fcalc(1) fcalc(end)]./1e9);
    ylim([-30 5]);
    xlabel('Frequency (GHz)');
    ylabel('Output Power (dB)');
    title('Output Power vs Frequency');
    set(gca,'FontSize',12,'FontWeight','bold');

    % --- Subplot 2: Output Power vs Position ---
    % The frequency used for this plot depends on the override setting.
    subplot(2, 1, 2);
    hold all;
    
    if ~isempty(single_frequency_override)
        % Use the specified override frequency for the position plot
        override_freq = single_frequency_override(1);
        [~, freq_idx_for_plot] = min(abs(fcalc - override_freq));
        plot_freq_ghz = fcalc(freq_idx_for_plot)/1e9;
        plot_title = sprintf('Output Power vs Position at Specified Freq (%.2f GHz)', plot_freq_ghz);
        main_title = sprintf('Full Sweep Analysis (Position Plot @ %.2f GHz, Ratio: %.3f)', plot_freq_ghz, current_ratio);
        save_filename = sprintf('SpecifiedFreq_Analysis_Ratio_%.3f.png', current_ratio);
    else
        % Find and use the optimal frequency for the position plot
        mode_to_optimize_idx = (harmonic_to_optimize_for + 1) / 2; % e.g., 3p is the 2nd mode
        output_power_of_harmonic = g_data(:, end, mode_to_optimize_idx);
        [~, freq_idx_for_plot] = max(output_power_of_harmonic);
        plot_freq_ghz = fcalc(freq_idx_for_plot) / 1e9;
        plot_title = sprintf('Output Power vs Position at Optimal Pump Freq for %dp (%.2f GHz)', harmonic_to_optimize_for, plot_freq_ghz);
        main_title = sprintf('Full Sweep Analysis with Optimization (Ratio: %.3f)', current_ratio);
        save_filename = sprintf('OptimalFreq_Analysis_Ratio_%.3f.png', current_ratio);
    end
    
    % Now, draw the position plot using the determined frequency index
    for h_idx = 1:ceil(maxHarmonic/2)
        plot(zcalc./twpa.len, smooth(g_data(freq_idx_for_plot,:,h_idx),1), 'Linewidth', 2);
    end
    legend(plotLegend);
    grid on;
    xlabel('Position (Normalized)');
    ylabel('Output Power (dB)');
    title(plot_title);
    set(gca,'FontSize',12,'FontWeight','bold');

    % Add a main title for the entire figure and save it
    sgtitle(f_handle, main_title);
    saveas(f_handle, fullfile(output_dir, save_filename));
    close(f_handle);
end

toc