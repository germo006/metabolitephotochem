function Cen = centroid2(concs, t)
%CENTROID computes the geometric centroids of the decay curves in the
%y-direction. This uses an inverted version of the decay equation and an
%analytically solved integral that, well, you'll sort of just need to
%believe.  
%   INPUTS
%       concs: n x 1 vector of concentrations
%       t:     n x 1 vector of durational times, where the first element is
%           treated as time zero. Literally, treated as a zero.
%   OUTPUTS
%       Cen:   n-1 x 1 vector of centroids

Co = concs(1);
Ct = concs(2:end);
Ct2 = Ct.^2;
t = t(2:end);

kp = -log(Co./Ct)./t;

Cen = Ct + 0.5.*(((Co^2)/2) - (Ct2./2) - (Ct2.*log(Co)) + (Ct2.*log(Ct)))./(...
    0.5*Co - kp.*t.*0.5.*Co);


end

