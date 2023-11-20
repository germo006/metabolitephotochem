function Wabs = WcylCompound(I, alph, r, eps) %, l)C,
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
%       Wabs: Light absorbed, mol photons L^-1 s^-1
% Because I'm tired and uninspired, this numerical integration will just
% use a trapezoidal sum.

%% Conversions
% nM to M
%C = C./1e9;
% m^2 to cm^2 
I = I./(100^2);
% m to cm
alph = alph./100;
r = r./10;
eps = eps./100;

prefactor = I.*eps./(alph);
x = 0:0.01*r:0.99*r;
toIntegrate = (1-10.^(-2.*alph*sqrt(r^2-x.^2)))./sqrt(r^2-x.^2);
integral = trapz(x, toIntegrate, 2);
Wabs = 1e3.*prefactor.*integral;

end

