% write velocity model to surf96 format file
% 
% created by Yang Zha 09/28/2012
% modified 10/17/2012 to include density in model
% modified 6/12/2023 to include Qp and Qs in model
% 
% model is nlayer*4 matrix [H Vp Vs Density]
% 		model(:,1): H(KM)
%		model(:,2): VP(KM/S)
%		model(:,3): VS(KM/S)
%		model(:,4): DENSITY(g/cm^3)	
% ---------- optional ----------
%       model(:,5): QP (optional)
%       model(:,6): QS (optional)
%
function fid=writemod_surf96(model,filename,varargin)
%
%
if nargin==2
    FREFP= 1; % default, reference frequency
    FREFS= 1; % default, reference frequency
    ETAP= 0; % default, eta
    ETAS= 0; % default, eta
elseif nargin==3
    FREFP = varargin{1};
    FREFS = varargin{1};
    ETAP= 0; % default, eta
    ETAS= 0; % default, eta
elseif nargin==4
    FREFP = varargin{1};
    FREFS = varargin{1};
    ETAP= varargin{2};
    ETAS= varargin{2}; 
end
%
% set up constant parameters
%rho=3.0; % constant density
%rho_water = 1.0; % water density
if size(model,2) == 4
    QP=1000*ones(size(model,1));
    QS=200*ones(size(model,1));
elseif size(model,2) == 6
    QP = model(:,5);
    QS = model(:,6);
else 
    error('Model must have either 4 or 6 columns');
end
fid = fopen(filename,'w');

	fprintf(fid,'MODEL.01\n');
	fprintf(fid,'model\n');
	fprintf(fid,'ISOTROPIC\n');
	fprintf(fid,'KGS\n');
%     fprintf(fid,'FLAT EARTH\n');
	fprintf(fid,'SPHERICAL EARTH\n');
	fprintf(fid,'1-D\n');
	fprintf(fid,'CONSTANT VELOCITY\n');
%     fprintf(fid,'VARIABLE VELOCITY\n');
	fprintf(fid,'LINE08\n');
	fprintf(fid,'LINE09\n');
	fprintf(fid,'LINE10\n');
	fprintf(fid,'LINE11\n');
	fprintf(fid,'      H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS    \n');
	for i=1:size(model,1)
		fprintf(fid,'     %7.4f',model(i,1)); % H(z)
		fprintf(fid,'     %7.4f',model(i,2));	% Vp
		fprintf(fid,'     %7.4f',model(i,3));	% Vs
		fprintf(fid,'     %7.4f',model(i,4));	% density
		fprintf(fid,'  	  %7.4f',QP(i)); % QP	
		if(model(i,3)<=0.01)
			fprintf(fid,'  	  %7.4f',0); % QS
		else
			fprintf(fid,'  	  %7.4f',QS(i)); % QS
        end
        % Q = Q ( f / f_refP )^eta
		fprintf(fid,'     %7.4f',ETAP); % ETAP
		fprintf(fid,'     %7.4f',ETAS); % ETAS
		fprintf(fid,'     %7.4f',FREFP); % FREFP
		fprintf(fid,'     %7.4f    ',FREFS); % FREFS
		fprintf(fid,'\n');
	end

	fclose(fid);
end
