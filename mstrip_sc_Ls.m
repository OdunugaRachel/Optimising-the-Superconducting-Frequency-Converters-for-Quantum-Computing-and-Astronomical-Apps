function [Ls, Cs, Z0, epeff, tandel, g2] = mstrip_sc_Ls( h, w, t1, t2, esub, ...
    eup, freq, Qsub, Qup, lambda1, lambda2)
% simple extension of mstrip to include extra inductance of superconducting 
% transmission line.  Have not put in frequency dependent effects, so only
% good when f << Delta/h.  Also have not included surface resistance.
%
%    h - substrate thickness (microns)
%    w - width of microstrip (microns)
%    t1 - conductor thickness - strip (microns)
%    t2 - conductor thickness - ground (microns)
%    esub - dielectric constant of substrate
%    eup - dielectric constant superstrate (usually = 1)
%    freq - frequency (GHz)
%    Qsub - dielectric quality factor of substrate (1/tan(delta))
%       delta is defined as the phase angle of the complex
%       dielectric constant of the substrate
%    Qup - dielectric quality factor of superstrate
%    lambda1 = magnetic penetration depth of strip
%    lambda2 = magnetic penetration depth of groundplane
%
%
%    Output parameters:
%
%    Z0 - characteristic impedance (ohms)
%    epeff - effective dielectric constant
%         i.e. wavelength = free space wavelength / sqrt(epeff)
%    tandel - effective loss tangent of line. Defined so that
%         beta = (2*Pi/wavelength)*(1 + i*tandel)
%         where beta is the complex propagation factor; i = sqrt(-1).
%    g2 - geometrical factor for loss calculations; = 1 for wide lines

   % call the normal microstrip function 
   [Zms, epeff, tandel, g2] = mstrip(h, w, t1, t2, esub, eup, freq, Qsub, Qup) ;

   % Convert to equivalent series and shunt elements
   cLight = 3e8;                    %m per second
   vph = cLight / sqrt( epeff);     %phase velocity on transmission line
   Ls = Zms / vph;                   %series inductance per unit length  H/m
   Cs = 1 / (Zms * vph);            %capacitance per unit length
   
   mu0 = pi*4e-7*1e-6;   %H/um
   Lsurf_cond = mu0 * lambda1 * coth( t1/ lambda1); %surface inductance of conductor
   Lsurf_gp =  mu0 * lambda2 * coth( t2/ lambda2);     %surface incutance of ground plane
   
   %find an effective width corresponding to current flow in the ground
   %plane - take this to be the pearl length for lambda > t2
   w_eff = max( [w lambda2.^2 / t2]);
   
   Ls = Ls + 1e6*Lsurf_cond/w + 1e6*Lsurf_gp / w_eff;            %corrected inductance per unit length
                                          %1e6 to covert to H/m - change
                                          %this back to  g2*Lsurf/w for
                                          %small lambda
   Z0 = sqrt( Ls/Cs);
   vph = 1/sqrt(Ls*Cs);
   epeff = (cLight / vph)^2; 
   
