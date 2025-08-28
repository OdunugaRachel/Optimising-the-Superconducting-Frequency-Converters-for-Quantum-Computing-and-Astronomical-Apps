function [Z0, epeff, tandel, g2] = mstrip( h, w, t1, t2, esub, eup, freq, Qsub, Qup)

%    Based on equations given by Hammerstad and Jensen, IEEE Microwave
%    Symposium proceedings, pp.407-409, 1980.
%
%    Lengths are converted to be in microns and frequencies in GHz.
%
%    It uses standard microwave engineering formulas (Hammerstad and Jensen) and 
%    includes corrections for conductor thickness, dispersion, and dielectric losses. 
%    The function is used to model how a microstrip line will behave in real circuits, especially at high frequencies.
%
%    Input parameters:
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


% Convert substrate dielectric constant to be relative to superstrate
%
      esub = esub/eup ;

      eta0 = 377; %ZVacuum/Ohm ;
      eta0 = eta0 / sqrt(eup) ;
      u = w/h ;

%     Calculate thickness corrections

      tnorm = t1/h ;                         % ?? need more thinking */
      ex = exp(sqrt(6.517*u)) ;
      coth = (ex+1./ex)./(ex-1./ex) ;
      delu1 = (tnorm/pi)*log(1. + 4*exp(1.)./(tnorm*coth.*coth)) ;
      ex = exp(sqrt(esub-1.)) ;
      cosh = 0.5*(ex+1./ex) ;
      delur = 0.5*(1. + 1./cosh).*delu1 ;
      u1 = u + delu1 ;
      ur = u + delur ;

%     Calculate characteristic impedance

      u = ur ;
      fu = 6. + (2.*pi-6.)*exp( -((30.666./u).^0.7528) ) ;
      Z01 = eta0.*log(fu./u + sqrt(1.+((2./u).*(2./u))))/(2.*pi) ;
      x1 = u.*u.*u.*u + (u/52.).*(u/52.) ;
      x2 = u.*u.*u.*u + 0.432 ;
      a = 1.+ (1./49.)*log(x1./x2) ;
      a = a + (1./18.7)*log(1.+(u/18.1).*(u/18.1).*(u/18.1)) ;
      b = 0.564*( ((esub-0.9)./(esub+3.)).^0.053 )  ;
      ee0 = (esub+1.)/2. + ((esub-1.)/2.).*( (1.+10./u).^(-a.*b) ) ;

      Z0 = Z01./sqrt(ee0) ;

%     calculate finite thickness corrections to effective diel. const.

      u = u1 ;
      fu = 6. + (2.*pi-6.)*exp( -1*( (30.666./u).^0.7528) ) ;
      Z02 = eta0.*log(fu./u + sqrt(1.+((2./u).*(2./u))))/(2.*pi) ;
      ee1 = ee0.*((Z02./Z01).*(Z02./Z01)) ;

%     corrections for dispersion

      G = (pi*pi/12.)*((esub-1.)./ee1).*sqrt(2.*pi*Z0./eta0) ;

%     Cutoff frequency in GHz of first TE mode ; here h in microns.

      fp = 397.887*Z0./h ;

%     dispersion-corrected dielectric constant

      epeff = esub - (esub-ee1)./(1.+G.*((freq./fp).*(freq./fp))) ;

%     dispersion-corrected characteristic impedance

      Z0 = Z0.*sqrt(ee1./epeff).*(epeff-1.)./(ee1-1.) ;

%     Attenuation - dielectric Losses

      q = (epeff - 1.)./(esub - 1.) ;
      Qd = ((1.-q) + q.*esub)./((1.-q)./Qup + q.*esub./Qsub) ;

      u = w./h ;
      Kfact = exp(-1.2.*((Z01./eta0).^0.7)) ;     %/* ?? need more thinking */
      g2 = Kfact;       %was g2 = 2*Kfact / w
                        %this seems a more sensible definition - P.Day

%     Assume no resistive losses - these will be taken into
%     account later

      Qc = 1e12; %MSTRIP_H_VERYBIG ;   /* just big ! */

%     Total attenuation

      Q0 = 1./(1./Qd + 1./Qc) ;
      tandel = 0.5./Q0 ;

      epeff = epeff .* eup ;     % we calculated relative to superstrate */
      
      
      