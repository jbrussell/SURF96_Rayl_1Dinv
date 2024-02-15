function [kernel, dispersion] = calc_anelastic_kernel96(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% calculate surface wave sensitivity kernel using srfker96 command from CPS  package :
%
% [kernel, dispersion] = calc_kernel96(model, period, wavetype, ifNorm, ifplot);
% model: surf96 format model: [thickness, vp, vs, rho, Qp, Qs];
% period : in seconds
% wavetype: 'L', 'R' or 'J' (joint);
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
%
% STEPS
% 1. make parameter file 
% 2. make model file and dispersion file 
% 3. use surf96  and srfker96 to calculate kernel
% 4. read in kernel using read_anelastic_kernel96.m
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
% jbrussell 6/12/2023: Significant rewrite to include anelastic kernels. This version
% reads in a new structure called "dispersion" that includes attenuation coefficient
% gamma as well as phv and grv with the influence of physical dispersion.
%

% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

if size(model,2)<6
    error('Must include Qp and Qs in model to generate anelastic kernels');
end

[nlayer,temp]=size(model);
Nper = length(vec_T);
% dc/da      dc/db      dU/da      dU/db      dc/dh      dU/dh     dc/dQai    dc/dQbi    dU/dQai    dU/dQbi    dg/dQai    dg/dQbi
kernel.dcdvs=nan*ones(nlayer,Nper);
kernel.dcdvp=nan*ones(nlayer,Nper);
kernel.dudvs=nan*ones(nlayer,Nper);
kernel.dudvp=nan*ones(nlayer,Nper);
kernel.dcdrho=nan*ones(nlayer,Nper);
kernel.dudrho=nan*ones(nlayer,Nper);
kernel.dcdQs=nan*ones(nlayer,Nper);
kernel.dcdQp=nan*ones(nlayer,Nper);
kernel.dudQs=nan*ones(nlayer,Nper);
kernel.dudQp=nan*ones(nlayer,Nper);
kernel.dgdQs=nan*ones(nlayer,Nper);
kernel.dgdQp=nan*ones(nlayer,Nper);
kernel.dgdQmu=nan*ones(nlayer,Nper);
kernel.dgdQkap=nan*ones(nlayer,Nper);
kernel.dQdQmu=nan*ones(nlayer,Nper);
kernel.dQdQkap=nan*ones(nlayer,Nper);

if nargin == 5
    nmode = 0; % default fundamental mode
    fref = 1; % default reference frequency
elseif nargin == 6
    nmode = varargin{1};
    fref = 1; % default reference frequency
elseif nargin == 7
    nmode = varargin{1};
    fref = varargin{2};
end

for ip = 1:Nper
    period = vec_T(ip);
    %% ================1. MAKE PARAMETER FILES====================
    make_par_surf96(wavetype,nmode); % make a control file for surf96
    % ================2. make model file and dispersion file ================
    writemod_surf96(model,'start.mod',fref); % write model to file

    fakedata = [period,3, 0.01]; % fake velocity entry for use by srfker96

    if(wavetype=='J') % if wavetype is joint L and R, write both to dispersion file, otherwise only write one
        error('Has not been updated to include Joint data')
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


    [kern_ip, disp_ip] = read_anelastic_kernel96('srfker96.txt',wavetype,nmode);
    if isempty(kern_ip)
        disp('Something went wrong calculating kernels');
        flds = fields(kern_ip);
        for ii = 1:length(flds)
            kernel.(flds{ii}) = [];
        end
        return
    end
    
    % Populate kernel structure for given frequency
    flds_kern = fields(kern_ip);
    for ii = 1:length(flds_kern)
        kernel.(flds_kern{ii})(:,ip) = kern_ip.(flds_kern{ii})(:);
    end
    
    % Populate dispersion structure for given frequency
    flds_disp = fields(disp_ip);
    for ii = 1:length(flds_disp)
        dispersion.(flds_disp{ii})(ip,1) = disp_ip.(flds_disp{ii});
    end
    
    % Choose kernel normalization
    if(ifnorm)
        for ii = 1:length(flds_kern)
            kernel.(flds_kern{ii})(:,ip) = kernel.(flds_kern{ii})(:,ip) ./ model(:,1);
        end
        t=['Depth sensitivity kernel',num2str(period),' S'];
    else
        t=['Layer sensitivity kernel',num2str(period),' S'];

    end
    % Set inf values to --> 0
    for ii = 1:length(flds_kern)
        Iinf = isinf(kernel.(flds_kern{ii})(:,ip));
        kernel.(flds_kern{ii})(Iinf,ip) = 0;
    end
    
    % Convert gamma(Qs,Qp) kernels to --> Q(Qs,Qp) kernels
    % gamma = omega / 2 / grv / Q
    omega = 2 * pi ./ dispersion.period(ip);
    kernel.dQdQs(:,ip) = (2 .* dispersion.grv(ip) ./ omega) .* kernel.dgdQs(:,ip);
    kernel.dQdQp(:,ip) = (2 .* dispersion.grv(ip) ./ omega) .* kernel.dgdQp(:,ip);
    
    % Convert gamma(Qs,Qp) kernels to --> gamma(Qmu, Qkappa) kernels
    % Derived following eq. 8 of Wang and Cai (2017) "A Method to Estimate Shear Quality Factor of Hard Rocks"
    % Qs^-1 = Qmu^-1
    % Qp^-1 = (1-4/3 * vs^2/vp^2)*Qkap^-1 + 4/3*(vs^2/vp^2)*Qmu^-1
    vp = model(:,2);
    vs = model(:,3);
    kernel.dgdQmu(:,ip) = kernel.dgdQs(:,ip) + kernel.dgdQp(:,ip).*(4/3 .* vs.^2 ./ vp.^2);
    kernel.dgdQkap(:,ip) = kernel.dgdQp(:,ip) .* (1 - 4/3 .* vs.^2 ./ vp.^2);
     
    % Convert gamma(Qmu, Qkappa) --> Q(Qmu, Qkappa) kernel
    omega = 2 * pi ./ dispersion.period(ip);
    kernel.dQdQmu(:,ip) = (2 .* dispersion.grv(ip) ./ omega) .* kernel.dgdQmu(:,ip);
    kernel.dQdQkap(:,ip) = (2 .* dispersion.grv(ip) ./ omega) .* kernel.dgdQkap(:,ip);

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
            plot(kernel.dcdvs(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVs');
            set(gca,'fontsize',20);
        elseif(wavetype=='R')
            subplot(1,4,1);
    %		plotlayermods(model(:,1),dcdvp,'-b'); 
            plot(kernel.dcdvp(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVp');
            set(gca,'fontsize',20);

            subplot(1,4,2);
    %		plotlayermods(model(:,1),dcdvs,'-b'); 
            plot(kernel.dcdvs(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('dC/dVs');
            set(gca,'fontsize',20);
            
            subplot(1,4,3);
    %		plotlayermods(model(:,1),dcdvp,'-b'); 
            plot(kernel.dQdQkap(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('d{Q^{-1}}/dQ^{-1}_{\kappa}');
            set(gca,'fontsize',20);

            subplot(1,4,4);
    %		plotlayermods(model(:,1),dcdvs,'-b'); 
            plot(kernel.dQdQmu(:,ip),z_mid,'-b'); axis ij

            title(t);
            xlabel('d{Q^{-1}}/dQ^{-1}_{\mu}');
            set(gca,'fontsize',20);
        end

    end
end

% Convert alpha (gamma) --> Q
omega = 2*pi./dispersion.period;
dispersion.qinv = 2 .* dispersion.grv .* dispersion.gamma ./ omega;

kernel.z = z_mid;

return





