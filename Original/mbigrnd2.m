function result = mbigrnd2( eps, delta, kbt, hbo)

%f1 = (exp(eps / kbt) + 1).^-1;
f2 = (exp((eps + hbo) / kbt) + 1).^-1;

result = (1 - 2*f2) .* abs(eps .* (eps + hbo) + delta.^2) ./ ...
      ( sqrt( delta^2 - eps.^2) .* sqrt( (eps + hbo).^2 - delta^2));


% function to calculate the frequency dependent response of a superconductor