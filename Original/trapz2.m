function result = trapz2( func, a, b, eps, npts, varargin)

	eps = (b - a) * eps;

	midp = 0.5 * (a+b);

	xpoints = a + logspace( log10(eps), log10( midp - a), npts / 2); 
	ypoints = feval( func, xpoints, varargin{:});

	i1 = trapint( ypoints, xpoints);
   
  	xpoints = b - logspace( log10(eps), log10( b - midp), npts / 2); 
	ypoints = feval( func, xpoints, varargin{:});
   
   i2 = -1*trapint( ypoints, xpoints);
   
   result = i1 + i2;

return


function trapi = trapint( y, x)

	trapi = 0.5 * sum( (x(2:end) - x(1:end-1)) .* (y( 2:end) + y( 1:end-1)));

return


% Numerical intergration method (Trapezoidal rule)
