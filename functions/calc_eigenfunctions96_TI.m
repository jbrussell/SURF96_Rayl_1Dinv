% Calculate surface wave eigenfunctions for a transversely isotropic (TI)
% model (TI analog of calc_eigenfunctions96.m)
%
%   eig = calc_eigenfunctions96_TI(vec_T,model,wavetype,nmode,fref)
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
% Output structure ([nlayer x nperiod] matrices; columns are NaN at periods
% where the requested mode does not exist):
%   eig.periods, eig.mode, eig.z (layer mid-point depth)
%   Rayleigh: eig.ur, eig.uz, eig.tr, eig.tz
%   Love:     eig.ut, eig.tt
%   In addition (not in the isotropic version):
%   eig.phv, eig.grv, eig.gamma, eig.amp, eig.ellip, and the phase-velocity layer
%   kernels eig.dcdvpv, eig.dcdvph, eig.dcdvsv, eig.dcdvsh, eig.dcdeta,
%   eig.dcdrho, eig.dcdh  (see read_der96_TI.m)
%
% Uses the CPS programs tprep96, tdisp96, tregn96/tlegn96 (via
% calc_der96_TI.m). The binary eigenfunction file is read directly
% (read_der96_TI.m) rather than the text output of tdpder96, which has only
% 3 significant figures. Spurious tdisp96 roots are rejected before the
% mode branch is selected (see calc_der96_TI.m).
%
% jbrussell 11/2026
%
function eig = calc_eigenfunctions96_TI(vec_T,model,wavetype,varargin)

if nargin==3
    nmode = 0; % default, fundamental mode
    fref = 1; % default, reference freq
elseif nargin==4
	nmode = varargin{1};
    fref = 1; % default, reference freq
elseif nargin==5
    nmode = varargin{1};
    fref = varargin{2};
end
wavetype = upper(wavetype);

% run tprep96, tdisp96 and tregn96/tlegn96 -DER, read the binary
% eigenfunction file and reject spurious roots (see calc_der96_TI.m)
eig = calc_der96_TI(vec_T,model,wavetype,nmode,fref);

% Get midpoint of layer
vec_h=model(:,1);
z1=cumsum(vec_h);z1=z1(:);
z2=[0;z1(1:end-1)];
z_mid=(z1+z2)/2;
eig.z = z_mid;

end
