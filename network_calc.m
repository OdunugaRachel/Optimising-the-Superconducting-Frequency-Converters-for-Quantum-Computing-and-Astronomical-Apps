clear 
%close all

tic

%g = geometry
%ms = microstrip parameters 
%finger = finger parameters
%twpa = travelling wave parametric amplifier

g.finger.l = 20;  %length of a capacitor 'finger', um add a bit to represent fringe effects
g.finger.w = 0.25; % dON't %.25;  %width of the microstrip line (and fingers)

% i thought that the width of the fingers didn't have to match the width of the microstrip line

g.finger.p = 0.87;    %period of the fingers (spacing between fingers), μm

% g.finger.l, g.finger.w and g.finger.p all directly influence the characterstic imedpedance of the microstrip

%% Don't change these parameters

g.ms.w = g.finger.w; %
g.ms.h = 0.15; %dielectric thickness
g.ms.eps = 10.3;  %dielectric epsilon
g.ms.eps_upper = 10.4;
g.ms.t1 = 0.030;  %microstrip conductor metal *thickness*
g.ms.t2 = 0.200;    %ground plane *thickness*

%% modulation parameters

%g.finger.modperiod = 50; %period of the modulation, μm
g.finger.modperiod = 50;
%   the length over+ which the finger pattern repeats.
g.finger.modamp = 10.8;% %modulation amplitude, μm, 
% 7 seems to be the max value that wokes without the frequency cutting off (with the intial gemotery parameters)
% controls how much the finger length varies.

unit_cells = 900;

L = g.finger.modperiod*unit_cells; % Total length of the device (μm)
% whay are we multiplying by 500 here? Number of unit cells
% the change in the length, make little differnece on the dispersion but has a big effect on the gain

filename = '0714TWPaX.mat';  
%loop over frequencies
% Array of frequencies (Hz) over which the simulation is performed.
f = linspace( .01, 50.01, 80000)*1e9; 

% 50 GHz is the upper limit of the simulation, this is the frequency at which the gain starts to drop off

%resistivities to consider
rhon = 480.*1e-8; %film resistivity, ohm.m
sn = 1 / rhon;  %normal state conductivity, (1/Ω·m), inverse of resistivity


%%

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

% Why is there an extra factor of a half in the gap calculation?
% What formula are we using to calculate the gap above? 
% (http://www.hyperphysics.phy-astr.gsu.edu/hbase/Solids/bcs.html)
% The gap is calculated as half the BCS energy gap at zero temperature.


lambda = zeros( 1, numel(f));  % London penetration depth (m), frequency-dependent.
% creating the array to store the London penertaration depth values

% Lperm: Per-unit-length inductance (H/m), frequency-dependent.
% Cperm: Per-unit-length capacitance (F/m), frequency-dependent.
% Z0: Characteristic impedance (Ω), frequency-dependent.

Lperm = lambda;  Cperm = lambda; Z0 = lambda; % I understand, we simply intialize these arrays

% z = 1:1:round(L/g.finger.modperiod/2);  % position array (m), half the period of the modulation.
% zmax = z(end); % last postition

% z and zmax are not used in the code

%% Calculate or Load Resistivity-based Z,vph

%calculate array of lambdas

for ii = 1:numel(f)

    s2 = mb2(gap, kbt, hbo(ii)) * sn;       %Mattis-Bardeen parameter, using mb2.m (function)
    lambda(ii) = 1 / sqrt(mu0 * s2 * omega(ii));  % calculating the london penertration depth for each frequency

    [Lperm(ii), Cperm(ii), Z0(ii), ~, ~, ~] = mstrip_sc_Ls( ms.h, ms.w, ms.t1, ms.t2, ...
        ms.eps, ms.eps_upper, 10, 1e9, 1e9, lambda(ii)*1e6, lambda(ii)*1e6);

        % mstrip_sc_Ls is a function that calculates the per-unit-length inductance, capacitance, and characteristic impedance of a superconducting microstrip line.
end

vph = 1 ./ sqrt(Lperm .* Cperm);  %m/s, phase velocity
IcOverIstar = 0.15;     % Ratio of critical current ot the star current
vph = vph./sqrt(1+IcOverIstar^2+IcOverIstar^4);

% I don't know where this is coming from

Dfinger = g.finger.p * 1e-6; %m, distance between fingers, period of the modulation in microns
Lfinger = g.finger.l * 1e-6; %m, Length of each finger (note:  2 fingers per section), in microns
Zref = 50;  %  Reference impedance for S-parameter calculation (Ω).
S21 = zeros(size(f)); % creating the S21 array

%%

% Retrieving the S-parameters for the device

ncell = round(L/g.finger.modperiod); % number of cells, based on the length and the modulation period

for ii = 1:numel(f)  % cycling through frequencies

    abcd_tot = [1 0; 0 1];  % Total ABCD matrix for the whole device. (idenity matrix)
    abcd = [1 0; 0 1];  % ABCD matrix for a single unit cell.
     
    n_unit_cell = g.finger.modperiod / g.finger.p ; % Number of unit cells per modulation period.
    
    %calculate abdc matrix for a unit cell
    for jj = 1:n_unit_cell 
        
        beta = 2*pi*f(ii) / vph(ii);
        betaf = 2*pi*f(ii) / vph(ii);

        % this is wave number
        
        Lfinger = (g.finger.w/2 + g.finger.l + g.finger.modamp * cos( 2*pi*(jj-0.5)*g.finger.p / g.finger.modperiod))* 1e-6; % length of each finger, with respect to the modulation period

        Zin_finger = -1i * Z0(ii) * cot( betaf * Lfinger); %  input impedence of the finger section

        %abcd matrix for shunt admittance - factor of 2 accounts for 2 fingers
        %per section
        abcd_finger = [1 0; 2/Zin_finger 1];
        % following theory of the refelection from a curved mirror

        %half length section - 
        %modify the length by a factor to try to represent the widening of
        %the line as it connects to the finger
        xfactor = (g.finger.p - g.finger.w)/g.finger.p;
        abcd_trl = [cos( beta * Dfinger * xfactor / 2)          1i * Z0(ii) * sin( beta * Dfinger * xfactor / 2); ...
            1i * sin( beta * Dfinger * xfactor / 2) / Z0(ii)    cos( beta * Dfinger * xfactor / 2)];

        %cascade the finger and TRL section S-matrices
        abcd = abcd * abcd_trl;
        abcd = abcd * abcd_finger;
        abcd = abcd * abcd_trl; 
    end

    %Combine into total abcd
    for n = 1:ncell
        abcd_tot = abcd_tot * abcd;
    end

    %Store new Sparamps for the position
    S = abcd2s(abcd_tot, Zref);
    S21(ii) = S(2, 1);
end

%%

%figure
%plot(f, angle(S21))

%%

dphi_2 = diff(unwrap(angle(S21)));
dphi = diff(angle(S21));
max_jump = max(abs(dphi));
max_jump_2 = max(abs(dphi_2));
disp(['Max phase jump (rad): ', num2str(max_jump)]);
disp(['Max phase jump - unwrapped (rad): ', num2str(max_jump_2)]);

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

%% Checking Dispersion Condition



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
% plot(f./1e9,-unwrap(angle(S21z(end,:))))
% plot(f,z2dB((S21z(end,:))),'Linewidth',2)
% xlim([0 30])
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
disp('________________________');
%%

L_in_cm = L./1e4;

save(filename, 'kperm', 'f', 'S21', 'len_meters', 'vph','Zfin')
% save('Target.mat', 'kperm', 'f', 'S21', 'len_meters', 'vph','Zfin')
save('Variables_Parameters.mat', 'Zfin', 'vph', 'g', 'L_in_cm')
saveas(f1, 'S21.png')
saveas(f2, 'dispersion.png')

toc