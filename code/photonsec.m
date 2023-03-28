function E = photonsec(I,l)
%PHOTONSEC converts a wavelength-specific irradiance (W m^-2) value and its
%corresponding wavelength (nm) to umol photons m^-2 s^-1

E = I.*l.*0.00836;

end

