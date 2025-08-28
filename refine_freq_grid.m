function [f_refined, S21, Lperm, Cperm, Z0, vph] = refine_freq_grid(f_init, g, ms, gap, kbt, mu0, sn, L, threshold)
% refine_freq_grid: Returns a refined frequency grid so that phase jump between points is < threshold (rad)
%
%   f_init: Initial frequency grid (Hz)
%   g: Geometry parameters struct
%   ms: Microstrip parameters struct
%   gap: Superconducting energy gap parameter (J/eV)
%   kbt: Thermal energy (J/eV)
%   mu0: Permeability of free space (H/m)
%   sn: Normal state conductivity (1/Ω·m)
%   L: Total length of the device (μm)
%   threshold: Maximum allowed phase jump (radians)
%
%   f_refined: Refined frequency grid (Hz)
%   S21: Transmission coefficient for the refined grid
%   Lperm: Per-unit-length inductance for the refined grid (H/m)
%   Cperm: Per-unit-length capacitance for the refined grid (F/m)
%   Z0: Characteristic impedance for the refined grid (Ω)
%   vph: Phase velocity for the refined grid (m/s)

    max_iter = 10;
    ev = 1.6e-19;

    for iter = 1:max_iter
        omega = 2*pi*f_init;
        hbo = 1e-34 * omega / ev;
        
        % Preallocate arrays for parfor
        Lperm = zeros(1, numel(f_init)); 
        Cperm = zeros(1, numel(f_init)); 
        Z0 = zeros(1, numel(f_init));

        % Use parfor for parallel execution if Parallel Computing Toolbox is available
        parfor ii = 1:numel(f_init)
            % 'current_lambda' is a temporary variable, local to each parallel worker.
            s2 = mb2(gap, kbt, hbo(ii)) * sn;
            current_lambda = 1 / sqrt(mu0 * s2 * omega(ii));
            [Lperm(ii), Cperm(ii), Z0(ii), ~, ~, ~] = mstrip_sc_Ls(ms.h, ms.w, ms.t1, ms.t2, ms.eps, ms.eps_upper, 10, 1e9, 1e9, current_lambda*1e6, current_lambda*1e6);
        end

        vph = 1 ./ sqrt(Lperm .* Cperm);
        IcOverIstar = 0.15;
        vph = vph./sqrt(1+IcOverIstar^2+IcOverIstar^4);

        Dfinger = g.finger.p * 1e-6;
        Zref = 50;
        S21 = zeros(size(f_init));

        ncell = round(L/g.finger.modperiod);
        
        % Use parfor for the main S21 calculation loop
        parfor ii = 1:numel(f_init)
            abcd_tot = [1 0; 0 1];
            abcd = [1 0; 0 1]; % ABCD matrix for a single unit cell

            n_unit_cell = g.finger.modperiod / g.finger.p;

            for jj = 1:n_unit_cell
                beta = 2*pi*f_init(ii) / vph(ii);
                betaf = 2*pi*f_init(ii) / vph(ii);
                Lfinger = (g.finger.w/2 + g.finger.l + g.finger.modamp * cos(2*pi*(jj-0.5)*g.finger.p / g.finger.modperiod))* 1e-6;
                Zin_finger = -1i * Z0(ii) * cot(betaf * Lfinger);
                abcd_finger = [1 0; 2/Zin_finger 1];
                xfactor = (g.finger.p - g.finger.w)/g.finger.p;
                abcd_trl = [cos(beta * Dfinger * xfactor / 2)          1i * Z0(ii) * sin(beta * Dfinger * xfactor / 2); ...
                    1i * sin(beta * Dfinger * xfactor / 2) / Z0(ii)    cos(beta * Dfinger * xfactor / 2)];
                abcd = abcd * abcd_trl;
                abcd = abcd * abcd_finger;
                abcd = abcd * abcd_trl;
            end

            for n = 1:ncell
                abcd_tot = abcd_tot * abcd;
            end

            %Store new Sparamps for the position
            S = abcd2s(abcd_tot, Zref);
            S21(ii) = S(2, 1);
        end

        dphi = diff(unwrap(angle(S21)));
        max_jump = max(abs(dphi));
        disp(['Max phase jump (rad): ', num2str(max_jump)]);

        idx = find(abs(dphi) >= threshold);
        if isempty(idx)
            f_refined = f_init;
            return; % Grid is sufficiently refined, assign output and exit the function.
        end

        N = 10; % Number of points to insert between each interval
        f_starts = f_init(idx).';
        f_ends = f_init(idx+1).';
        steps = (f_ends - f_starts) / (N + 1);
        new_points_matrix = f_starts + (1:N) .* steps;
        new_points = new_points_matrix(:).';
        
        all_points = [f_init, new_points];
        tol = 1; % Hz
        all_points = sort(all_points);
        keep = [true, diff(all_points) > tol];
        f_init = all_points(keep);
    end
    warning('Adaptive grid refinement reached max iterations. Some phase jumps may still exceed threshold.');
    f_refined = f_init; % Ensure f_refined is assigned if max_iter is reached
end