clear 
%close all

tic

%g = geometry
%ms = microstrip parameters 
%finger = finger parameters
%twpa = travelling wave parametric amplifier

g.finger.l = 26;  %length of a capacitor 'finger', um add a bit to represent fringe effects
g.finger.w = 0.25; % dON't %.25;  %width of the microstrip line (and fingers)

% i thought that the width of the fingers didn't have to match the width of the microstrip line

g.finger.p = 2;    %period of the fingers (spacing between fingers), μm

% g.finger.l, g.finger.w and g.finger.p all directly influence the characterstic imedpedance of the microstrip

%% Don't change these parameters

g.ms.w = g.finger.w; %
g.ms.h = 0.19; %dielectric thickness
g.ms.eps = 10.3;  %dielectric epsilon
g.ms.eps_upper = 10.4;
g.ms.t1 = 0.035;  %microstrip conductor metal *thickness*
g.ms.t2 = 0.200;    %ground plane *thickness*

%% modulation parameters

%g.finger.modperiod = 50; %period of the modulation, μm
g.finger.modperiod = 110;
%   the length over+ which the finger pattern repeats.
g.finger.modamp = 4; %modulation amplitude, μm, 
% 7 seems to be the max value that wokes without the frequency cutting off (with the intial gemotery parameters)
% controls how much the finger length varies.

unit_cells = 500;

L = g.finger.modperiod*unit_cells; % Total length of the device (μm)
% the change in the length, make little differnece on the dispersion but has a big effect on the gain

filename = '0714TWPaX.mat';  
%loop over frequencies
% Array of frequencies (Hz) over which the simulation is performed.
f = linspace( .01, 50.01, 60000)*1e9; 

% 50 GHz is the upper limit of the simulation, this is the frequency at which the gain starts to drop off

%resistivities to consider
rhon = 260.*1e-8; %film resistivity, ohm.m
sn = 1 / rhon;  %normal state conductivity, (1/Ω·m), inverse of resistivity

%loop over position
ms = g.ms;

ev = 1.6e-19;  %J/ev
mu0 = pi*4e-7; %H/m
omega = 2*pi*f;
hbo = 1e-34 * omega / ev;

% odd way of defining the Bolzmann constant and temperature

kb = 1.38e-23;  %J/K
T = 0.01;  %operating temperature, K
kbt = kb * T / ev; % If temperture is changing, I would made this a function of T

Tc = 14.5;    % Critical temperature of the superconductor (K)
gap = 0.5 * 3.5 * Tc * kb / ev;   % Superconducting energy gap parameter (J/eV).


% Adaptive frequencygrid refinements

phase_jump_threshold = 0.85*pi; % threshold fo the phase jumps
[f, S21, Lperm, Cperm, Z0, ~] = refine_freq_grid(f, g, ms, gap, kbt, mu0, sn, L, phase_jump_threshold);


omega = 2*pi*f;

%%

f_GHz = f / 1e9;

uth = -1*unwrap(angle(S21(:)));
for fj = 2:length(f)
    while uth(fj)<uth(fj-1)
        uth(fj:end) = uth(fj:end) + 2*pi;
    end
end

slope = (uth(2)-uth(1))./(f(2)-f(1));
intercept = uth(1) - slope*f(1);
uth = uth - intercept;

len_meters = L./1e6;
kperm = uth/len_meters;

p = polyfit(f(f<0.1e9), kperm(f<0.1e9), 1); % 
vph = 2*pi/p(1);


%% Plot Results

f1 = figure(1);
hold on
plot(f./1e9,abs(S21))
grid on
grid minor
title('s21 vs Frequency')
xlabel('Frequency (GHz)')
ylabel('S_{21}')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1000 100 1500 1000])


f2 = figure(2);
hold on
plot(f./1e9,kperm,'Linewidth',2)

grid on
grid minor
title('Dispersion vs Frequency')
xlabel('Frequency (GHz)')
ylabel('Dispersion (m^{-1})')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
set(gcf,'Position',[1000 100 1500 1000])

Ctot = Cperm.*(2*g.finger.l + g.finger.p)./g.finger.p;

Zfin = real(1./(vph*Ctot)); 
 
disp(['Impedance (Z) = ',num2str(Zfin(1)),' Ohms'])
disp(['Phase Velocity (vph) = ',num2str(vph./3e8),' c'])
disp(['Length (L) = ',num2str(L./1e4),' cm'])
disp(['Frequency Array Size = ',num2str(size(f))])
disp('________________________');
%%

L_in_cm = L./1e4;

save(filename, 'kperm', 'f', 'S21', 'len_meters', 'vph','Zfin')
save('Variables_Parameters.mat', 'Zfin', 'vph', 'g', 'L_in_cm', 'unit_cells')
saveas(f1, 'S21.png')
saveas(f2, 'dispersion.png')

toc