
clear

%% Initialize TWPA

% load paramp_NbTiN_JPL_mK_short879_259_281_x0.875_taper_highf.mat
load C:\Users\klimovich\Documents\MATLAB\TWPA\0714TWPaX.mat

% data = importdata("C:\Users\klimovich\Downloads\E08_06_750Cells_CAngDeg.csv");

twpa = createTWPA;

% twpa.fsim = data(:,1).*1e9;
% twpa.ksim = -data(:,2)/360*2*pi/twpaLEN;

% data = importdata("C:\Users\klimovich\Downloads\E08_06_750Cells_S21.csv");

% twpa.gsim = -data(:,3)/twpaLEN;
% twpa.gsim = twpa.gsim - twpa.gsim;


twpa.fsim = f;
twpa.ksim = kperm;
twpa.gsim = -log(abs(transpose(S21)));

% Sanitize NaNs
twpa.ksim(isnan(twpa.ksim)) = max(twpa.ksim);
twpa.gsim(isnan(twpa.gsim)) = -100;

twpa.pumpF = 6.3e9;
% twpa.pumpF2 = 0e9;
twpa.Istar = 4.5*1000;
twpa.Ip = twpa.Istar*0.10;
twpa.Idc = twpa.Istar*0.0;
    
twpa.len = 110e-6*879*1;
twpa.betanl = 1;

%% Pick Modes

% modes = {'p','3p','5p','7p','9p'};
modes = {'p','s','i','p+s','p+i'};

for i=1:length(modes)
    twpa = twpa.addMode(modes{i},4);
end

% Alternative for dual-pump: with each mode being a row of 
% [(# of pump1 photons) (# of signal photons) (# of pump 2 photons)]
% so the example below is "p1, s, p2, 2p1 - s, 2p2 - s
% twpa.modes = [1 0 0; 0 1 0; 0 0 1; 2 -1 0; 0 -1 2];

twpa.I0 = zeros(length(twpa.modes),1);
twpa.I0(1) = twpa.Ip;
twpa.I0(2) = twpa.Ip.*1e-3;

%% Compute Gain

% Frequencies over which to calculate gain
fcalc = 0.01e9:0.05e9:19.01e9;
fcalc = removeHarmonics(fcalc,[twpa.pumpF 3*twpa.pumpF],twpa.Idc);

% Positions to sample
zcalc = 0:0.001:twpa.len;

g = zeros(length(fcalc), length(zcalc), length(twpa.modes));
Iend = zeros(length(fcalc), length(zcalc), length(twpa.modes));

for ii = 1:length(fcalc)
    wn = twpa.modes(:,1)*twpa.pumpF + twpa.modes(:,2)*fcalc(ii);
%    If using dual pump instead use:
%    wn = twpa.modes(:,1)*twpa.pumpF + twpa.modes(:,2)*fcalc(ii) + twpa.modes(:,3)*twpa.pumpF2;
    S21 = exp((-twpa.g(wn.') + 1i.*twpa.k(wn.')).*twpa.len);

    Y = solveCME(fcalc(ii),zcalc,twpa);

    g(ii,:,:) = 20*log10(abs(Y(:,:).*S21./twpa.I0(2)));  
    disp(ii/length(fcalc))
end

%% Plot Results

figure(1)
cmap = colormap('turbo');
hold all
% for i=1:length(twpa.modes)
plot(fcalc./1e9,smooth(g(:,end,2),1),'Linewidth',2)
% plot(fcalc./1e9,smooth(g(:,end,4),1),'Linewidth',2)
% plot(zcalc,g(1,:,1),'Linewidth',2)
% plot(zcalc,g(1,:,2),'Linewidth',2)
% plot(zcalc,g(1,:,3),'Linewidth',2)
% plot(zcalc,g(1,:,4),'Linewidth',2)
% plot(zcalc,g(1,:,5),'Linewidth',2)
% end
grid on
% xlim([fcalc(1) fcalc(end)]./1e9)
ylim([0 30])
xlabel('Frequency (GHz)')
ylabel('Gain (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
drawnow

% figure(2)
% cmap = colormap('turbo');
% hold all
% plot(smooth(g(1,:,1),1),'Linewidth',2)
% plot(smooth(g(1,:,4),1),'Linewidth',2)
% plot(smooth(g(1,:,5),1),'Linewidth',2)
% plot(smooth(g(1,:,6),1),'Linewidth',2)
% plot(smooth(g(1,:,7),1),'Linewidth',2)
% legend('p','3p','5p','7p','9p')
% % plot(smooth(g(1,:,8),1),'Linewidth',2)
% % plot(smooth(g(1,:,9),1),'Linewidth',2)
% grid on
% % xlim([fcalc(1) fcalc(end)]./1e9)
% % ylim([-20 0])
% xlabel('Position')
% ylabel('Gain (dB)')
% set(gca,'FontSize',16)
% set(gca,'FontWeight','bold')
% set(gcf,'Position',[1500 100 1500 1000])
% drawnow



