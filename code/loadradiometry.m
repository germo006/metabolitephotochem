function [l_SunTest,W,E] = loadradiometry
% LOADRADIOMETRY.m: this is just a function to upload a file of radiometry
% (that is, wavelength-dependent irradiance values) which were measured by
% Danielle Freeman in the SunTest Atlas XLS+. 
%   The file to be loaded is "Radiometry_RulerSunTest_forMATLAB.xlsx". This
%   is the only option, so this function requires no input. 
data = readtable("../datasets/Radiometry_RulerSunTest_forMATLAB.xlsx");
l_SunTest = data.l;
% The data are spatial, so I will create a data cube called W that contains
% the wavelenth as the long axis and the x,y dimensions as the primary
% dimensions. 
W = zeros(3,4,height(data));
W(1,1,:) = data.A1; W(1,2,:) = data.A2; W(1,3,:) = data.A3; W(1,4,:) = data.A4;
W(2,1,:) = data.B1; W(2,2,:) = data.B2; W(2,3,:) = data.B3; W(2,4,:) = data.B4;
W(3,1,:) = data.C1; W(3,2,:) = data.C2; W(3,3,:) = data.C3; W(3,4,:) = data.C4;

% Now a version of this that contains the same values, converted to photon
% flux in uE m^-2
E = zeros(3,4,height(data));
converter = @(E) photonsec(E,l_SunTest);
E(1,1,:) = converter(data.A1); E(1,2,:) = converter(data.A2);
E(1,3,:) = converter(data.A3); E(1,4,:) = converter(data.A4);
E(2,1,:) = converter(data.B1); E(2,2,:) = converter(data.B2);
E(2,3,:) = converter(data.B3); E(2,4,:) = converter(data.B4);
E(3,1,:) = converter(data.C1); E(3,2,:) = converter(data.C2);
E(3,3,:) = converter(data.C3); E(3,4,:) = converter(data.C4);
end

