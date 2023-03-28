function Wabs = Wcyl(I, alph, r) %, l)
%WABS converts incoming irradiance (resolved to 1 nm) and screening factors
%obtained from spectrophotometry to absorbed radiation at each wavelength
%for a cylinder on its side with radius r and length l
%   INPUT
%       I:    Irradiance, umol photons m^-2 s^-1,   n x m array
%       alph: screening factor, m^-1,               n x 1 array
%       r:    Cynlinder radius, mm,                 single value
%       l:    Cylinder length, mm,                  single value UNUSED
%
%   OUTPUT
%       Wabs: Light absorbed, umol photons m^-3 s^-1
% Because I'm tired and uninspired, this numerical integration will just
% use a trapezoidal sum.
prefactor = I.*10^9./(pi*r^2);
x = 0:0.1:r;
toIntegrate = 10.^(-2.*alph*sqrt(r^2-x.^2));
integral = trapz(x, toIntegrate, 2);
Wabs = prefactor.*(2*r - 2.*integral);

end

