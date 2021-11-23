function [outmodel, phv_fwd]= invdispR(vec_T,phv,phvstd,grv,grvstd,startmodel,mohodepth,waterdepth,seddepth,niteration,lambda,ddlay)
% INVERT 1D velocity model from Rayleigh wave dispersion, using SURF96
% INVERT FOR SHEAR VELOCITY AND KEEP THICKNESS FIXED
% ALSO KEEP VP/VS FIXED, DENSITY FROM VP 
%:INPUT: 
%	vec_T: vector of periods
%	phv: phase velocity data
%	phvstd : error of phase velocity
%	grv: group velocity data
%	grvstd : error of group velocity
%	startmodel : nlayer x 4 array of data:
%		model(:,1): H(KM) layer thicknesses
%		model(:,2): VP(KM/S)
%		model(:,3): VS(KM/S)
%		model(:,4): DENSITY(g/cm^3)
%	mohodepth: Depth to moho, used to break smoothing constraint
%	waterdepth: depth to seafloor from surface (used to determine where to allow sharp change)
%   seddepth: depth to sediments, used ot break smoothing constraint
%   lambda: controls damping in the inversion (larger values = more damped)
%   ddlay: controls smoothing across discontinuities
%
%OUTPUT: 
%	outmodel: inveted 1D model nlayer x 3 array of data:
%	phv_fwd: predicted phase velocity dispersion of final inverted model.
%
%	modified by Ge Jin, jinwar@gmail.com
%	2013-02-21
%
%Created by Yang Zha 09/27/2012
% modified 10/17/2012 to include density in model
%
% jbrussell 5/27/2020: Modify slightly how discontinuities are handled. Currently
% set up to handle 3 discontinuities (water-sediments, sediment-crust,
% crust-mantle). Added damping and smoothing parameters as inputs... found
% that default values were very unstable (underdamped).
%
% jbrussell 6/1/2020: Add timeout string to surf96 inversion call to avoid
% hanging. Get gtimeout on mac by installing coreutils with homebrew
% >$ brew install coreutils

% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

%%(0) generate a sobs.d file in the current directory (constant inversion parameters)
%make_par_surf96('R'); % make a control file for surf96
display('===========1D Rayleigh Dispersion Inversion============');
% clean up

system('surf96 39'); 
system('rm start.mod');
system('rm disp_obs.dsp');
%% (1) write starting model to file in current directory  

display('Writing model file..');

writemod_surf96(startmodel,'start.mod');

%% (2) write dispersion data to file
%
%
display('Writing dispersion data...');

if length(phv)==length(vec_T)

	Cdata(:,1)=vec_T(:);
	Cdata(:,2)=phv(:);
	Cdata(:,3)=phvstd(:);

	writedisp_surf96(Cdata,'disp_obs.dsp'); %write phase velocity data
end
if length(grv)==length(vec_T)
	Udata(:,1)=vec_T(:);
	Udata(:,2)=grv(:);
	Udata(:,3)=grvstd(:);

	writedisp_surf96(Udata,'disp_obs.dsp','R','U',1) %write group velocity data
end


%% (3) set smoothing constraints, allow sharp change of velocity at moho and ocean surface


vec_z=cumsum(startmodel(:,1)); % convert from layer thickness to depth

nl=length(startmodel(:,1)); % number of layers

% mohodepth=waterdepth+h_crust; % moho depth


% Smoothing scheme: permit large smoothing at seafloor and moho. larger change in the crust than in the mantle
%
%
display('Setting smoothing profile...');

system('surf96 36 1'); % set smoothign type

if waterdepth > 0
%     ind_sf = find(abs(vec_z-waterdepth)==min(abs(vec_z-waterdepth))); % find layer number correspinding to seafloor
    ind_h2o = find(startmodel(:,2)==1.5);
    if ~isempty(ind_h2o)
        ind_sf = ind_h2o(end);
    	cmdtemp = ['surf96 31 ',num2str(ind_sf),' ',num2str(ddlay)]; % large change at seafloor
    	system(cmdtemp);
    % 	cmdtemp = ['surf96 48 1 0']; % turn off the smoothing at that boundary
    % 	system(cmdtemp);
        cmdtemp = ['surf96 48 ',num2str(ind_sf),' 0']; % turn off the smoothing at that boundary
        system(cmdtemp);
    end
end

if seddepth>0
	ind_sed = find(abs(vec_z-seddepth)==min(abs(vec_z-seddepth))); % find layer number correspinding to moho depth
	cmdtemp = ['surf96 31 ',num2str(ind_sed(1)),' ',num2str(ddlay)]; % large change at seafloor
    system(cmdtemp);
    cmdtemp = ['surf96 48 ',num2str(ind_sed(1)),' 0']; % turn off the smoothing at that boundary
	system(cmdtemp);
end

if mohodepth>0
	ind_moho = find(abs(vec_z-mohodepth)==min(abs(vec_z-mohodepth))); % find layer number correspinding to moho depth
	cmdtemp = ['surf96 31 ',num2str(ind_moho(1)),' ',num2str(ddlay)]; % large change at seafloor
    system(cmdtemp);
    cmdtemp = ['surf96 48 ',num2str(ind_moho(1)),' 0']; % turn off the smoothing at that boundary
	system(cmdtemp);
end

% % Test what happens if everything is considered a "boundary"
% for ii = 1:length(vec_z)
%     cmdtemp = ['surf96 31 ',num2str(ii),' ',num2str(ddlay)]; % large change at seafloor
%     system(cmdtemp);
%     cmdtemp = ['surf96 48 ',num2str(ii),' 0']; % turn off the smoothing at that boundary
% 	system(cmdtemp);
% end

%% Adjust how Vp (and therefore density) is handled
% If vp is held fixed, so too will density. If vp/vs is held fixed, density
% will vary according to vp, but it's not clear to me exactly how...
is_fixVp = 0; % 1: fix vp / 0: fix vp/vs
if is_fixVp
    for ilay = 1:nl
        system(['surf96 30 ',num2str(ilay),' 0']);
    end
    
%     % Fix Vp in the mantle only
%     ind_mantle = find(vec_z > mohodepth);
%     for ilay = ind_mantle(:)'
%         system(['surf96 30 ',num2str(ilay),' 0']);
%     end
end


%% (4) Perform inversion
%
%####
%	start the first inversion with a slightly higher damping
%	do avoid an overshoot in the first model estimate
%####
display('Start inversion....');
display(['do ', num2str(niteration), ' iterations ..']);

system(['surf96 32 ',num2str(lambda*10),' > log_surf96.txt']); % define damping lambda
system([timeoutstr,'surf96 37 2 1 2 6 >> log_surf96.txt']); % start with 2 iterations
%####
%	do 10 inversions
%####
display('Inverting ..');

% Do inversion
cmdtemp=[timeoutstr,'surf96 32 ',num2str(lambda),' 37 ',num2str(niteration),' 1 2 6 >> log_surf96.txt'];
system(cmdtemp);

% cmdtemp=['surf96 32 1 37 ',num2str(niteration),' 1 2 6 >> log_surf96.txt'];
% system(cmdtemp);

display('Writing log file..');

system('surf96 45 >> log_surf96.txt'); % write smoothing info into log

system('surf96 17 >> log_surf96.txt'); % write dispersion of final model info into log

system([timeoutstr,'surf96 28 final.mod']); % output model

system([timeoutstr,'surf96 1 27 temp.dsp']);
%

system('surf96 39'); % clean up


%%(5) Read output model from file
display('Getting inverted model from:  final.mod..');

outmodel = readmod_surf96('final.mod');
% read dispersion from tmep file
data=readdisp_surf96('temp.dsp');

if(~isempty(data));
	phv_fwd=data(:,2);
else
	phv_fwd=[]
end


display('=================Inversion Finished========================');

return
end



