% Calculate Rayleigh wave admittance following Ruan et al. (2014) for a
% transversely isotropic (TI) model (TI analog of calc_admittance96_Ruan14.m)
%
% 	model: TI model [H Vpv Vsv rho Vph Vsh eta] or [H Vpv Vsv rho Vph Vsh eta Qp Qs]
%          (see writemod_surf96_TI.m; use convert_iso2TI.m to convert an
%          isotropic model). The water layer(s) must have Vsv = 0.
% 	vec_T: period
%   nmode: mode branch of interest (0 = fund.)
%
% jbrussell 11/2026
%
function adm = calc_admittance96_Ruan14_TI(vec_T,model,nmode)

% Calculate vertical displacement and stress eigenfunctions
eig = calc_eigenfunctions96_TI(vec_T,model,'R',nmode);

% Extract eigenfunctions at seafloor and calculate admittance
ind_h2o = find(model(:,3)==0);
iseafloor = ind_h2o(end)+1;
admittance = eig.uz(iseafloor,:) ./ eig.tz(iseafloor,:); % km/GPa
admittance = admittance * 1e-6; % convert km/GPa --> m/Pa

adm.admittance = admittance;
adm.periods = eig.periods;

end
