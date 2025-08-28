function Is = solveCME(fcalc, zcalc, twpa)

    % Create array of frequencies to consider
    wn = zeros(1,length(twpa.modes));
    for i=1:length(twpa.modes)
        wn(i) = twpa.modes(i,1)*twpa.pumpF + twpa.modes(i,2)*fcalc;
    end
    wn_orig = wn;
        
    % Remove duplicates and unphysical terms
    wn = wn(wn > 0);
    wn = unique(wn, 'stable');

    % Set initial conditions
    Y0 = zeros(1,numel(wn));
    
    for i=1:length(wn)
        duplicates = find(wn(i) == wn_orig);
        Y0(i) = twpa.I0(duplicates(1));
    end

    beta = twpa.k(wn.');
    gamma = twpa.g(wn.');

    % Solve CME
    pamp = @(x,y) CMEgen(wn, y, beta, gamma, x, twpa.Istar, twpa.Idc);
    [pos, As] = ode45(pamp,[0 twpa.len],Y0,odeset(RelTol=1e-6));

    Is = zeros(length(zcalc),length(wn_orig));
    for i=1:length(wn)
        % Restore duplicates
        duplicates = find(wn(i) == wn_orig);
        for j=1:length(duplicates)
            Is(:,duplicates(j)) = interp1(pos,As(:,i),zcalc);
        end
    end
%     Is = zeros(length(zcalc),length(wn_orig));

    % Recombine duplicates
%     for i=1:length(wn_orig)
%         ind = find(wn == wn_orig(i));
%         Is(:,i) = As(ind);
% 
%         
%         
%     end


end