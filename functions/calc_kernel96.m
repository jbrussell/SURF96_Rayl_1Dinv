function [dcdvs, dcdvp, dudvs, dudvp, z_mid, dcdrho, dudrho] = calc_kernel96(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% calculate surface wave sensitivity kernel using srfker96 command from CPS  package :
%
% [dc/dvs dc/dvp du/dvs du/dvp] = calc_kernel96(model, period, wavetype, ifNorm, ifplot);
% model: surf96 format model: [thickness, vp, vs, rho];
% period : in seconds
% wavetype: 'L', 'R' or 'J' (joint);
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
%
% STEPS
% 1. make parameter file 
% 2. make model file and dispersion file 
% 3. use surf96  and srfker96 to calculate kernel
% 4. read in kernel using read_kernel96.m
% 5. normalize and plot if specified
%
% Yang Zha 11/06/2013
%
% jbrussell 6/5/2020: Note, if using the kernels for velocity inversion,
% set ifnorm=0. Kernels are scaled such that input data 
% is in units of km/s and the output model is velocity perurbation in
% units of km/s (not percent!).
%
% jbrussell 6/1/2020: Add timeout string to surf96 inversion call to avoid
% hanging. Get gtimeout on mac by installing coreutils with homebrew
% >$ brew install coreutils
%
% jbrussell 11/20/2022: Include density (rho) kernels in addition to Vs and Vp.
% Also update to allow for overtone kernels to be output
% nmode=0 for fund. mode; nmode=1 for 1st overtone, etc...
%

% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

[nlayer,temp]=size(model);
Nper = length(vec_T);
dcdvs=nan*ones(nlayer,Nper);
dcdvp=nan*ones(nlayer,Nper);
dudvs=nan*ones(nlayer,Nper);
dudvp=nan*ones(nlayer,Nper);
dcdrho=nan*ones(nlayer,Nper);
dudrho=nan*ones(nlayer,Nper);

if nargin == 5
    nmode = 0; % default fundamental mode
elseif nargin == 6
    nmode = varargin{1};
end

for ip = 1:Nper
    period = vec_T(ip);
    %% ================1. MAKE PARAMETER FILES====================
    make_par_surf96(wavetype,nmode); % make a control file for surf96
    % ================2. make model file and dispersion file ================
    writemod_surf96(model,'start.mod'); % write model to file

    fakedata = [period,3, 0.01]; % fake velocity entry for use by srfker96

    if(wavetype=='J'); % if wavetype is joint L and R, write both to dispersion file, otherwise only write one

        writedisp_surf96(fakedata,'disp_obs.dsp','L','U',1,nmode);
        writedisp_surf96(fakedata,'disp_obs.dsp','R','U',0,nmode);
    else
        overwrite=1;
        writedisp_surf96(fakedata,'disp_obs.dsp',wavetype,'U',overwrite,nmode);
    %     writedisp_surf96(fakedata,'disp_obs.dsp');
    end

    %================ 3. use surf96  and srfker96 to calculate kernel% ================================
    system(['surf96 39']);
    [~, ~] = system([timeoutstr,'surf96 1']);
    system([timeoutstr,'srfker96 > srfker96.txt']);
    system('surf96 39');
    system('rm start.mod');
    system('rm disp_obs.dsp');
    %%================ READ KERNEL FILE AND ARRANGE============================


    kernel= read_kernel96('srfker96.txt',wavetype,nmode);
    if isempty(kernel)
        disp('Something went wrong calculating kernels');
        dcdvs = [];
        dcdvp = [];
        dudvs = [];
        dudvp = [];
        dcdrho = [];
        dudrho = [];
        z_mid = [];
        return
    end
    dcdvs(:,ip)=kernel(:,1);
    dcdvp(:,ip)=kernel(:,2);
    dudvs(:,ip)=kernel(:,3);
    dudvp(:,ip)=kernel(:,4); 
    dcdrho(:,ip)=kernel(:,5); 
    dudrho(:,ip)=kernel(:,6); 

    if(ifnorm)
        dcdvs(:,ip) = dcdvs(:,ip)./model(:,1);
        dcdvp(:,ip) = dcdvp(:,ip)./model(:,1);
        dudvs(:,ip) = dudvs(:,ip)./model(:,1);
        dudvp(:,ip) = dudvp(:,ip)./model(:,1);
        dcdrho(:,ip) = dcdrho(:,ip)./model(:,1);
        dudrho(:,ip) = dudrho(:,ip)./model(:,1);
        t=['Depth sensitivity kernel',num2str(period),' S'];
    else
        t=['Layer sensitivity kernel',num2str(period),' S'];

    end
    dcdvs(isinf(dcdvs(:,ip)),ip) = 0;
    dcdvp(isinf(dcdvp(:,ip)),ip) = 0;
    dudvs(isinf(dudvs(:,ip)),ip) = 0;
    dudvp(isinf(dudvp(:,ip)),ip) = 0;
    dcdrho(isinf(dcdrho(:,ip)),ip) = 0;
    dudrho(isinf(dudrho(:,ip)),ip) = 0;

    % Get mid point of layer
    vec_h=model(:,1);
    z1=cumsum(vec_h);z1=z1(:);
    z2=[0;z1(1:end-1)];
    z_mid=(z1+z2)/2;


    if(ifplot)
        % plotting sensitivity 
        figure;hold on
        if(wavetype=='L')
    %		plotlayermods(model(:,1),dcdvs,'-b'); 
            plot(dcdvs(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVs');
            set(gca,'fontsize',20);
        elseif(wavetype=='R');
            subplot(1,2,1);
    %		plotlayermods(model(:,1),dcdvp,'-b'); 
            plot(dcdvp(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVp');
            set(gca,'fontsize',20);

            subplot(1,2,2);
    %		plotlayermods(model(:,1),dcdvs,'-b'); 
            plot(dcdvs(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVs');
            set(gca,'fontsize',20);
        end

    end
end

return





