
clear

%% Initialize TWPA

% load paramp_NbTiN_JPL_mK_short879_259_281_x0.875_taper_highf.mat
load 0714TWPaX.mat

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

% twpa.pumpF = 6.3e9;
% twpa.pumpF2 = 0e9;
twpa.Istar = 4.5*1000;
twpa.Ip = twpa.Istar*0.08;
twpa.Idc = twpa.Istar*0.0;
    
twpa.len = 110e-6*879*1;
twpa.betanl = 1;

data = importdata("0714TWPaX.mat")
%% Pick Modes

maxHarmonic = 9;

twpa.modes = [1 0;   % fundamental
               3 0]; % third harmonic

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
fcalc = 0.1e9:0.1e9:5.1e9;

% Positions to sample
zcalc = 0:0.001:twpa.len;

g = zeros(length(fcalc), length(zcalc), length(twpa.modes));
Iend = zeros(length(fcalc), length(zcalc), length(twpa.modes));

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

pumpGrid  = linspace(0.5e9,2e9,50);    % e.g. 0.5→2 GHz
modAmpGrid = linspace(2,6,25);         % e.g. 2 μm → 6 μm

best = struct('ratio',-inf,'f',0,'modamp',0);
for modAmp = modAmpGrid
  % rebuild your unit-cell ABCD once here…
  for fp = pumpGrid
    twpa.pumpF = fp;
    % re-run your solveCME block for this single fp
    % [gout] = your existing code that ends up in g(end,:,1:2)
    G1 = squeeze( g(end,:,1) );  % 1st harmonic gain along z
    G3 = squeeze( g(end,:,2) );  % 3rd harmonic
    r  = max(G3) - max(G1);     % “how much hotter is 3rd vs 1st?”
    if r>best.ratio
      best.ratio  = r;
      best.f      = fp;
      best.modamp = modAmp;
      best.G1     = max(G1);
      best.G3     = max(G3);
    end
  end
end

fprintf('Best f_p = %.2f GHz, modAmp = %.2f μm → G3=%.1fdB, G1=%.1fdB\n',...
        best.f/1e9,best.modamp,best.G3,best.G1);
%% Plot Results
% close all

figure(1)
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
xlabel('Position')
ylabel('Power (dB)')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1500 100 1500 1000])
drawnow

index = length(fcalc);

figure(2)
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
drawnow

