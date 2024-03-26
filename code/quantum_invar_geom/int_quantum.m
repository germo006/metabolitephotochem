function Y = int_quantum(C,I,eps,alph,Q,t,lambda,r)
% Computes an integral on absorbed wavelengths lambda. 
%   INPUTS
%       I:      nx1 vector of irradiance, umol m-2 s-1 nm-1
%       eps:    nx1 molar extinction coefficients, M-1 m-1
%       Q:      1x1 quantum yield
%       kp:     1x1, estimated first-order decay rate, h-1
%       tf:     1x1, final time in h
%       lambda: nx1 range of integration wavelengths, nm
%       r:      radius of the tube, mm
%
%   OUTPUT
%       Y: an estimate of the value C0-Ct
%
%

Cbar = centroid2(C,t); % Find the geometric mean acting concentrations on 
                       % all time periods.

Wabs = WcylCompound(I, alph, r, Cbar, eps); % Integrate radiation absorbed
    % by the analyte over the time period and tube dimensions by
    % wavelength.

Y = Q.*trapz(lambda', Wabs, 1); % Sum over wavelengths.