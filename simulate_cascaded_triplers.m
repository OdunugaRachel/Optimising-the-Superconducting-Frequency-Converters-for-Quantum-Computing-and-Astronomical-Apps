clear
tic

%% STAGE 1: First Tripler Simulation

filename1 = '0714TWPaX_First.mat'; % Data file for the first tripler
data1 = load(filename1);

twpa1 = createTWPA;
maxHarmonic1 = 3;
twpa1.modes = [1 0];
for i=3:2:maxHarmonic1
    twpa1.modes = cat(1, twpa1.modes, [i 0]);
end

twpa1.Istar = 4.5*1000;
twpa1.Ip = twpa1.Istar*0.13;
twpa1.Idc = twpa1.Istar*0.0;
twpa1.betanl = 1;
twpa1.I0 = zeros(length(twpa1.modes),1);
twpa1.I0(1) = twpa1.Ip;


twpa1.fsim = data1.f;
twpa1.ksim = data1.kperm;
twpa1.gsim = -log(abs(transpose(data1.S21)));
twpa1.len = (110e-6*879*1);

% Sanitize NaNs robustly
twpa1.ksim(isnan(twpa1.ksim)) = max(twpa1.ksim);
twpa1.gsim(isnan(twpa1.gsim)) = -100;



pumpF1 = 3.6e9; % Pump frequency for the first device
fcalc1 = 0.1e9:0.01e9:5.1e9; 
zcalc1 = 0:0.0001:twpa1.len;

twpa1.pumpF = pumpF1;

fprintf('Input Pump: %.2f GHz, Current: %.3f A\n', pumpF1/1e9, abs(twpa1.I0(1)));


% --- Run Simulation for Stage 1 ---

Y1 = solveCME(pumpF1, zcalc1, twpa1);

% --- Extract Output for Stage 2 ---
I_out_1p = Y1(end, 1);
I_out_3p = Y1(end, 2); % Generated 3p harmonic

fprintf('This 3p current will be the pump for Stage 2.\n\n');


%% STAGE 2: Second Tripler (pumped by the output of the first)

filename2 = '0714TWPaX_Second.mat'; % Data file for the second tripler

data2 = load(filename2);

twpa2 = createTWPA;
maxHarmonic2 = 3;
twpa2.modes = [1 0];
for i=3:2:maxHarmonic2
    twpa2.modes = cat(1, twpa2.modes, [i 0]);
end
pumpF2 = 3 * pumpF1;
Ip2 = I_out_3p; 
twpa2.I0 = zeros(length(twpa2.modes),1);
twpa2.I0(1) = Ip2;

twpa2.fsim = data2.f;
twpa2.ksim = data2.kperm;
twpa2.gsim = -log(abs(transpose(data2.S21)));
twpa2.len =(110e-6*879*1)/5; % Use length from the data file

% Sanitize NaNs robustly
twpa2.ksim(isnan(twpa2.ksim)) = max(twpa2.ksim);
twpa2.gsim(isnan(twpa2.gsim)) = 100;

twpa2.Istar = 4.5*1000;
twpa2.Idc = twpa2.Istar*0.0;
twpa2.betanl = 1;

twpa2.pumpF = pumpF2;

fprintf('Input Pump: %.2f GHz, Current: %.3f A (complex: %.3f + %.3fi)\n', pumpF2/1e9, abs(twpa2.I0(1)), real(twpa2.I0(1)), imag(twpa2.I0(1)));

% --- Run Simulation for Stage 2 ---
fcalc2 = 8e9:0.01e9:12.5e9;
zcalc2 = 0:0.0001:twpa2.len;
Y2 = solveCME(pumpF2, zcalc2, twpa2);


%% Plot Results

f12 = figure(12);
set(gcf, 'Position', [100, 100, 1200, 500]);

% Define a reference current for dB conversion (initial pump power)
I_ref = twpa1.Ip;

% --- Plot for Stage 1 ---
subplot(1, 2, 1);
hold on;
plot(zcalc1, 20*log10(abs(Y1(:,1)) / abs(I_ref)), 'LineWidth', 2, 'DisplayName', sprintf('Pump (%.2f GHz)', pumpF1/1e9));
plot(zcalc1, 20*log10(abs(Y1(:,2)) / abs(I_ref)), 'LineWidth', 2, 'DisplayName', sprintf('3rd Harmonic (%.2f GHz)', 3*pumpF1/1e9));
hold off;
grid on;
xlabel('Position along device (mm)');
ylabel('Power (dB)');
title('Stage 1: First Tripler');
legend('show', 'Location', 'best');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylim([-80, 5]); % Set y-axis limits for better visualization

% --- Plot for Stage 2 ---
subplot(1, 2, 2);
hold on;
plot(zcalc2/max(zcalc2), 20*log10(abs(Y2(:,1)) / abs(I_ref)), 'LineWidth', 2, 'DisplayName', sprintf('Pump (%.2f GHz)', pumpF2/1e9));
plot(zcalc2/max(zcalc2), 20*log10(abs(Y2(:,2)) / abs(I_ref)), 'LineWidth', 2, 'DisplayName', sprintf('9p Harmonic (%.2f GHz)', 3*pumpF2/1e9));
hold off;
grid on;
xlabel('Position along device (mm)');

ylabel('Power (dB)');
title('Stage 2: Second Tripler');
legend('show', 'Location', 'best');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylim([-80, 5]); % Set y-axis limits for better visualization

sgtitle('Cascaded Tripler Power Evolution', 'FontSize', 16, 'FontWeight', 'bold');
drawnow;

saveas(gcf, 'Cascaded_Tripler_Power_Evolution.png');

g1 = zeros(length(fcalc1), length(zcalc1), length(twpa1.modes));  % gain array, creating empty array to store the gain values
g2 = zeros(length(fcalc2), length(zcalc2), length(twpa2.modes));  % gain array, creating empty array to store the gain values

for ii = 1:length(fcalc1)
    wn = twpa1.modes(:,1)*fcalc1(ii);
    twpa1.pumpF = fcalc1(ii);
    S21_linear_prop = exp((-twpa1.g(wn.') + 1i.*twpa1.k(wn.')).*twpa1.len);

    Y = solveCME(fcalc1(ii),zcalc1,twpa1);

    g1(ii,:,:) = 20*log10(abs(Y(:,:).*S21_linear_prop./twpa1.I0(1)));  
    disp(ii/length(fcalc1))
end

for ii = 1:length(fcalc2)
    wn = twpa2.modes(:,1)*fcalc2(ii);
    twpa2.pumpF = fcalc2(ii);
    S21_linear_prop = exp((-twpa2.g(wn.') + 1i.*twpa2.k(wn.')).*twpa2.len);

    Y = solveCME(fcalc2(ii),zcalc2,twpa2);

    g2(ii,:,:) = 20*log10(abs(Y(:,:).*S21_linear_prop./twpa2.I0(1)));  
    disp(ii/length(fcalc2))
end

f14 = figure(14);
subplot(1, 2, 1);

hold on
plotLegend = {};
for i=1:ceil(maxHarmonic1/2)
    plot(fcalc1./1e9,smooth(g1(:,end,i),1),'Linewidth',2)
    plotLegend = [plotLegend, {[num2str(1 + 2*(i-1)),'p']}];
end
legend(plotLegend, 'Location', 'best')
grid on
xlim([fcalc1(1) fcalc1(end)]./1e9)
ylim([-30 5])
xlabel('Frequency (GHz) ')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
title('Tripler 1: Output Power vs Frequency')
drawnow

figure(14);
subplot(1, 2, 2);
hold on
plotLegend = {};
for i=1:ceil(maxHarmonic1/2)
    plot(fcalc2./1e9,smooth(g2(:,end,i),1),'Linewidth',2)
    plotLegend = [plotLegend, {[num2str(1 + 2*(i-1)),'p']}];
end
legend(plotLegend, 'Location', 'best')
grid on
xlim([fcalc2(1) fcalc2(end)]./1e9)
ylim([-30 5])
xlabel('Frequency (GHz) ')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
title('Tripler 2: Output Power vs Frequency')
drawnow

toc