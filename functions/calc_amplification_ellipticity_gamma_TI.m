% Calculate surface wave dispersion, amplification, ellipticity, and
% attenuation for a transversely isotropic (TI) model
% (TI analog of calc_amplification_ellipticity_gamma.m)
%
%   disper = calc_amplification_ellipticity_gamma_TI(vec_T,model,wavetype,nmode,fref)
%
% 	model: TI model [H Vpv Vsv rho Vph Vsh eta] or [H Vpv Vsv rho Vph Vsh eta Qp Qs]
%          (see writemod_surf96_TI.m; use convert_iso2TI.m to convert an
%          isotropic model). If Qp and Qs are not provided, Qp=1000 and Qs=200
%          are used, as in the isotropic version.
% 	vec_T: period
%   wavetype: Rayleigh (R) or Love (L)
%   nmode: mode branch of interest (0 = fund.) (OPTIONAL, default 0)
%   fref: reference frequency of Q model in Hz (OPTIONAL, default 1)
%
% Output structure (same fields as the isotropic version):
%   disper.period, disper.phv, disper.grv, disper.gamma
%   Rayleigh: disper.A_R (amplification), disper.RZ (ellipticity)
%   Love:     disper.A_L (amplification)
% Unlike the isotropic version, all requested periods are returned (in the
% order requested), with NaN where the mode does not exist.
%
% The isotropic version reads the ASCII table written by sdpegn96. Here the
% same quantities (C, U, GAMMA, ENERGY = AR/AL, ELLIPTICITY) are read from
% the binary tregn96/tlegn96 -DER output, which has full single precision
% and lets spurious tdisp96 roots be rejected (see calc_der96_TI.m). The
% values are identical to those written by "tdpegn96 -C -U -G -A0 -E -ASC".
%
% jbrussell 11/2026
%
function disper = calc_amplification_ellipticity_gamma_TI(vec_T,model,wavetype,varargin)

if nargin==3
    nmode = 0; % default, fundamental mode
    fref = 1; % default, reference frequency
elseif nargin==4
	nmode = varargin{1};
    fref = 1; % default, reference frequency
elseif nargin==5
    nmode = varargin{1};
    fref = varargin{2}; % reference frequency
end
wavetype = upper(wavetype);

eig = calc_der96_TI(vec_T,model,wavetype,nmode,fref);

disper.period = eig.periods;
disper.phv = eig.phv;     % phase velocity km/s
disper.grv = eig.grv;     % group velocity km/s
disper.gamma = eig.gamma; % attenuation coefficient (1/km)
if strcmp(wavetype,'R')
    disper.A_R = eig.amp;  % Rayleigh amplification
    disper.RZ = eig.ellip; % ellipticity
elseif strcmp(wavetype,'L')
    disper.A_L = eig.amp;  % Love amplification
end

end
