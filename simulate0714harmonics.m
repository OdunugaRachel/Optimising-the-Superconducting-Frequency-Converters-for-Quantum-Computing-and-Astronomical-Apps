 
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
    
twpa.len = 110e-6*879*1;   % length of the the device, I believe this is fixed
% where does this length come from and whta does each factor represent?
twpa.betanl = 1;   %nonlinearity parameter, affects the strength of the nonlinear effects 
                   % %Josesphon Junction = 0.5


%% Pick Modes

maxHarmonic = 3;  %max harmonic order to simulate
 
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
fcalc = 0.1e9:0.01e9:9e9;  % pump frequencies

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
% close all

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
title('Output Power vs Frequency')
drawnow

index = length(fcalc);

f4 = figure(4);
cmap = colormap('turbo');
hold all
for i=1:ceil(maxHarmonic/2)
    plot(zcalc./twpa.len,smooth(g(index,:,i),1),'Linewidth',2)
end
legend(plotLegend)
% plot(smooth(g(1,:,8),1),'Linewidth',2)
% plot(smooth(g(1,:,9),1),'Linewidth',2)
grid on
% xlim([fcalc(1) fcalc(end)]./1e9)
% ylim([-20 0])
xlabel('Position')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
title('Output Power vs Position')
drawnow

% %% Sum of squared amplitudes (energy conservation check)
% 
% sum_sq = zeros(length(fcalc),1);
% for ii = 1:length(fcalc)
%     amps = squeeze(g(ii,end,:)); % dB values
%     amps_lin = 10.^(amps/20);   % convert dB to linear amplitude
%     sum_sq(ii) = sum(amps_lin.^2); % sum of squares
% end
% 
% f_sum = figure(5);
% plot(fcalc./1e9, sum_sq, 'k-', 'LineWidth', 2)
% grid on
% xlabel('Frequency (GHz)')
% ylabel('Sum of squared amplitudes')
% title('Sum of squared amplitudes at output (should be constant if energy conserved)')
% set(gca,'FontSize',16)
% set(gca,'FontWeight','bold')
% set(gcf,'Position',[1500 100 1500 1000])
% drawnow
% saveas(f_sum, 'SumSquaredAmplitudes.png')

saveas(f3, 'Frequency_Power.png')
saveas(f4, 'Harmonics.png')

toc