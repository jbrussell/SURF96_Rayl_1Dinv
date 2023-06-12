% function of calculating rayleigh wave dispersion 
% 	model: [thickness vp vs]
% 	vec_T: period 
% 	data: [phv];
%   nmode: 0=fund. 1=1st, ... etc. (OPTIONAL, defaults to fund.)
%   disp_type: 'C'=phase veloc. 'U'=group veloc (OPTIONAL, defaults to C)
%
% jbrussell 11/20/2022: Update to allow calculation of overtone data
% 
function phv = dispR_surf96(vec_T,model,varargin)

% jbrussell 6/16/2020: Add timeout string to surf96 inversion call to avoid
% hanging. Get gtimeout on mac by installing coreutils with homebrew
% >$ brew install coreutils
% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

system('surf96 39'); % clean up
system('rm start.mod');
system('rm temp.dsp');

if nargin==2
    nmode = 0; % default, fundamental mode
    disp_type = 'C'; % default, phase velocity
elseif nargin==3
	nmode = varargin{1};
    disp_type = 'C'; % default, phase velocity
elseif nargin==4
    nmode = varargin{1};
    disp_type = upper(varargin{2});
    if ~(strcmp(disp_type,'C') || strcmp(disp_type,'U'))
        error('disp_type must be C or U');
    end
end

% make surf96 par file
make_par_surf96('R',nmode);

% Make dummy data
datatemp(:,1) = vec_T(:);
datatemp(:,2) = 3;
datatemp(:,3)  =0.1;

writedisp_surf96(datatemp,'disp_obs.dsp','R',disp_type,1,nmode); % write temp data into dispersion file 
writemod_surf96(model,'start.mod');

% run surf96 

[~, log] = system([timeoutstr,'surf96 1 27 temp.dsp']);
% read dispersion from tmep file
data=readdisp_surf96('temp.dsp');

if(~isempty(data))
	phv=data(:,2);
else
	phv=[];
end
%
system('surf96 39'); % clean up

end


