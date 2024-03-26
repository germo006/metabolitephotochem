function Wabs = WcylCompound(I, alph, r, Cbar, eps) %, l)
%WABS converts incoming irradiance (resolved to 1 nm) and screening factors
%obtained from spectrophotometry to absorbed radiation at each wavelength
%for a cylinder on its side with radius r and length l
%   INPUT
%       I:    Irradiance, umol photons m^-2 s^-1,   n x m array
%       alph: screening factor, m^-1,               n x 1 array
%       r:    Cynlinder radius, mm,                 single value
%       l:    Cylinder length, mm,                  single valueNOT USED
%       C:    Concentration, nM,                    1 x m array
%       eps:  Molar absorbance, m^-1 M^-1,          n x 1 array
%
%   OUTPUT
%       Wabs: Light absorbed, mol photons L^-1 s^-1
% Because I'm tired and uninspired, this numerical integration will just
% use a trapezoidal sum.

%% Conversions
% nM to M
Cbar = Cbar./1e9;
% m^2 to cm^2 
I = I./(100^2); %umol photons s^-1 cm^-2
% m to cm
alph = alph./100; %cm^-1
r = r./10;  %cm
eps = eps./100; % cm^-1 M^-1

prefactor = 4.*I.*(eps*Cbar')./(pi.*alph.*r^2);
x = 0:0.01*r:0.99*r;
toIntegrate = 1-10.^(-2.*alph*sqrt(r^2-x.^2));
integral = trapz(x, toIntegrate, 2);
Wabs = 1e3.*prefactor.*integral;

end

