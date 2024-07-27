% Calculate Rayleigh wave admittance following Ruan et al. (2014)
% 	model: [thickness vp vs rho]
% 	vec_T: period 
%   nmode: mode branch of interest (0 = fund.)
%
% jbrussell 7/2024
% 
function adm = calc_admittance96_Ruan14(vec_T,model,nmode)

% Calculate vertical displacement and stress eigenfunctions
eig = calc_eigenfunctions96(vec_T,model,'R',nmode);

% Extract eigenfunctions at seafloor and calculate admittance
ind_h2o = find(model(:,3)==0);
iseafloor = ind_h2o(end)+1;
admittance = eig.uz(iseafloor,:) ./ eig.tz(iseafloor,:); % km/GPa
admittance = admittance * 1e-6; % convert km/GPa --> m/Pa

adm.admittance = admittance;
adm.periods = eig.periods;

end


