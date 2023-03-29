function Y = int_quantum(W,eps,Q,kp,tf,lambda)
% Computes a double integral to time tf on wavelengths lambda. 
%   INPUTS
%       W:      nx1 vector of irradiance, umol m-3 s-1 nm-1
%       eps:    nx1 molar extinction coefficients, M-1 cm-1
%       Q:      nx1 quantum yields
%       kp:     1x1, estimated first-order decay rate, h-1
%       tf:     1x1, final time in h
%       lambda: nx1 range of integration wavelengths, nm
%
%   OUTPUT
%       Y: an estimate of the value C0-Ct/(2.3*C0)
%
%

% Unit conversions
eps_m_nM = eps.*100./1e9;  % Now nM-1 m-1
W_h_nM = W.*3600.*1000;     % Now nM m-3 h-1 nm-1

t = 0:1/60:tf;

yh = (W_h_nM.*eps_m_nM.*Q)*exp(-kp.*t);

Y = trapz(t, trapz(lambda', yh, 1),2);