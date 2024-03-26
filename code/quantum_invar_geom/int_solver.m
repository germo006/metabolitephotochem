function dydt = int_solver(W,eps,Q,lambda)
% Computes an integral on wavelengths lambda to get a degradation rate at
% a specific time. 
%   INPUTS
%       W:      nx1 vector of irradiance, umol m-2 s-1 nm-1
%       eps:    nx1 molar extinction coefficients, M-1 cm-1
%       Q:      nx1 quantum yields
%       kp:     1x1, estimated first-order decay rate, h-1
%       C:      1x1, concentration in nM
%       lambda: nx1 range of integration wavelengths, nm
%
%   OUTPUT
%       dydt
%
%
% Unit conversions
eps_m_nM = eps.*100./1e9;  % Now nM-1 m-1
W_h_nM = W.*3600.*1000;     % Now nM m-2 h-1 nm-1

yh = W_h_nM.*eps_m_nM.*Q;
dydt = 2.3.*trapz(lambda', yh);