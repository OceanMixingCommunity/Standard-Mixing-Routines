function Pt = inertial_model(epsi,k,veldir)
% Required inputs
% veldir=1,2,or 3 for longitudinal, transverse or vertical direction
% k=Spectral  wavenumber k (rad/m)

% Outputs:
%     Pt: theoretical spectra in (m/s)^2/(rad/m)
if veldir==1
     Aj=1;
else
     Aj=4/3;
end

alp=1.5*(18/55)*Aj;% kte for theoretical spectrum

%Pt= alp*f.^(-5/3)*(epsi.*U).^(2/3);% theoretical spectrum in rad/s.
Pt= alp*k.^(-5/3)*epsi.^(2/3);% theoretical spectrum in (m/s)^2/rad/m.

end