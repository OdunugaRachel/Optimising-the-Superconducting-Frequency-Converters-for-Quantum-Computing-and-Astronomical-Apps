function result = mb2(delta, kbt, hbomega)

result = zeros( length(kbt), 1);

for ii = 1:length(kbt)
result(ii) = trapz2( 'mbigrnd2', delta(ii) - hbomega, delta(ii), ...
   1e-10, 6000, delta(ii), kbt(ii), hbomega) / hbomega;
end


% Mattis-Bardeen parameter calculation, integrating the Mattis-bardeen formula
% delta is the superconducting energy gap (J)
% kbt is the thermal energy (J)