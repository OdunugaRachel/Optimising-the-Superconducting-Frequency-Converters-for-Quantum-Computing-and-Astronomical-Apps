function [dB] = z2dB(z)
    dB = 10*log10(real(z).^2 + imag(z).^2);
end