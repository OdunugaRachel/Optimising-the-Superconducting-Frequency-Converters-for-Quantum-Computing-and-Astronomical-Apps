function dI = CMEgen(wn, In, kn, gamma, z, Istar, Idc)

    N = length(wn);
    en = exp((-gamma + 1i*kn)*z);
    gn = -gamma + 1i*kn;

    % Create conjugate arrays of frequencies and currents
    ws = [-wn wn];
    Is = [conj(In.*en) In.*en];
    
    % Create tensors
    Tijk = Is(:)*Is(:).';
    Mijk = ws(:)+ws(:).';

    % Compute 3WM terms
    dI = zeros(N,1);
    if Idc ~= 0
        for x=1:N
            dI(x) = -2*kn(x)*kn(x)*Idc*sum((wn(x) == Mijk).*Tijk,"all")/(8*Istar^2*gn(x)*en(x));
        end
    end
    
    % Expand to 4WM tensors
    Tijk = reshape(Tijk(:).*Is(:).',[2*N,2*N,2*N]);
    Mijk = reshape(Mijk(:)+ws(:).',[2*N,2*N,2*N]);

    % Compute 4WM terms
    for x=1:N
        dI(x) = dI(x) - kn(x)*kn(x)*sum((wn(x) == Mijk).*Tijk,"all")/(24*Istar^2*gn(x)*en(x));
    end

end