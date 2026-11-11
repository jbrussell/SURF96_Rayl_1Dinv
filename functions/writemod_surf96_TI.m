% Write a transversely isotropic (TI) velocity model to a CPS MODEL96 file
% (TI analog of writemod_surf96.m)
%
% fid = writemod_surf96_TI(model, filename)
% fid = writemod_surf96_TI(model, filename, fref)
% fid = writemod_surf96_TI(model, filename, fref, eta_Q)
%
% model is nlayer x 7 (or nlayer x 9) matrix:
%       model(:,1): H(KM)
%       model(:,2): VPV(KM/S)  vertically travelling P velocity  sqrt(C/rho)
%       model(:,3): VSV(KM/S)  SV velocity                       sqrt(L/rho)
%       model(:,4): DENSITY(g/cm^3)
%       model(:,5): VPH(KM/S)  horizontally travelling P velocity sqrt(A/rho)
%       model(:,6): VSH(KM/S)  SH velocity                       sqrt(N/rho)
%       model(:,7): ETA        F/(A-2L)  (=1 for an isotropic layer)
% ---------- optional ----------
%       model(:,8): QP
%       model(:,9): QS
%
% Columns 1-4 have the same meaning as the isotropic model used by
% writemod_surf96.m ([H Vp Vs rho]), so an isotropic model is converted with
% model_TI = [model(:,1:4), model(:,2), model(:,3), ones(nlayer,1)]
% (see convert_iso2TI.m). As in writemod_surf96.m, if QP and QS are not
% provided, QP = 1000 and QS = 200 are written. A Q value of 0 is treated by
% CPS as an elastic (Q^-1 = 0) layer.
%
% The CPS TI model file stores VPF = sqrt(F/rho) rather than eta; here
% F = eta*(A - 2L), i.e. VPF = sqrt(eta*(VPH^2 - 2*VSV^2)).
%
% jbrussell 11/2026
%
function fid=writemod_surf96_TI(model,filename,varargin)
%
if nargin==2
    FREFP= 1; % default, reference frequency
    FREFS= 1; % default, reference frequency
    ETAP= 0; % default, Q frequency exponent
    ETAS= 0; % default, Q frequency exponent
elseif nargin==3
    FREFP = varargin{1};
    FREFS = varargin{1};
    ETAP= 0; % default
    ETAS= 0; % default
elseif nargin==4
    FREFP = varargin{1};
    FREFS = varargin{1};
    ETAP= varargin{2};
    ETAS= varargin{2};
end
%
nlay = size(model,1);
if size(model,2) == 7
    QP=1000*ones(nlay,1);
    QS=200*ones(nlay,1);
elseif size(model,2) == 9
    QP = model(:,8);
    QS = model(:,9);
else
    error('TI model must have either 7 or 9 columns: [H Vpv Vsv rho Vph Vsh eta (Qp Qs)]');
end

H   = model(:,1);
VPV = model(:,2);
VSV = model(:,3);
RHO = model(:,4);
VPH = model(:,5);
VSH = model(:,6);
ETA = model(:,7);

% VPF = sqrt(F/rho), with F = eta*(A-2L)
F_over_rho = ETA .* (VPH.^2 - 2*VSV.^2);
if any(F_over_rho < 0)
    error('writemod_surf96_TI: eta*(Vph^2-2*Vsv^2) < 0 in some layers; F<0 cannot be written to a CPS TI model (VPF is stored)');
end
VPF = sqrt(F_over_rho);

fid = fopen(filename,'w');

	fprintf(fid,'MODEL.01\n');
	fprintf(fid,'model\n');
	fprintf(fid,'TRANSVERSE ISOTROPIC\n');
	fprintf(fid,'KGS\n');
%     fprintf(fid,'FLAT EARTH\n');
	fprintf(fid,'SPHERICAL EARTH\n');
	fprintf(fid,'1-D\n');
	fprintf(fid,'CONSTANT VELOCITY\n');
	fprintf(fid,'LINE08\n');
	fprintf(fid,'LINE09\n');
	fprintf(fid,'LINE10\n');
	fprintf(fid,'LINE11\n');
	fprintf(fid,'      H(KM)  VPV(KM/S)  VSV(KM/S) RHO(GM/CC)         QP         QS       ETAP       ETAS      FREFP      FREFS\n');
	fprintf(fid,'             VPH(KM/S)  VSH(KM/S)  VPF(KM/S)\n');
	for i=1:nlay
		if(VSV(i)<=0.01)
			qs = 0; % fluid layer
		else
			qs = QS(i);
		end
		% line 1: H VPV VSV RHO QP QS ETAP ETAS FREFP FREFS
		fprintf(fid,'     %7.4f     %7.4f     %7.4f     %7.4f     %10.4f     %10.4f     %7.4f     %7.4f     %7.4f     %7.4f\n', ...
			H(i), VPV(i), VSV(i), RHO(i), QP(i), qs, ETAP, ETAS, FREFP, FREFS);
		% line 2: VPH VSH VPF
		fprintf(fid,'                 %7.4f     %7.4f     %9.6f\n', VPH(i), VSH(i), VPF(i));
	end

	fclose(fid);
end
