function [dcdvs, dcdvp, dudvs, dudvp, z_mid, dcdrho, dudrho, dcdh, dudh, phv, grv] = calc_kernel96(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% calculate surface wave sensitivity kernel using srfker96 command from CPS  package :
%
% [dcdvs, dcdvp, dudvs, dudvp, z_mid, dcdrho, dudrho, dcdh, dudh, phv, grv] = ...
%       calc_kernel96(model, period, wavetype, ifNorm, ifplot, nmode, fref);
% model: surf96 format model: [thickness, vp, vs, rho];
% period : in seconds
% wavetype: 'L', 'R' or 'J' (joint);
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
% nmode: 0 for fund. mode; 1 for 1st overtone, etc. (OPTIONAL, default 0)
% fref: reference frequency (Hz) of the Q model (OPTIONAL, default 1; only
%       used if Qp and Qs are given)
%
% OUTPUT
% dcdvs, dcdvp, dudvs, dudvp: Vs and Vp kernels of phase (c) and group (U)
%       velocity (from srfker96)
% z_mid: layer mid-point depth
% dcdrho, dudrho: density kernels (from sregn96/slegn96 -DER, see
%       calc_density_kernel96.m)
% dcdh, dudh: the layer-thickness ("dc/dh", "dU/dh") kernels output by
%       srfker96
% phv, grv: elastic phase and group velocity
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
% jbrussell 11/2026:
%  - dcdrho/dudrho previously held columns 7-8 (Rayleigh) or 5-6 (Love) of
%    the srfker96 output, which are dc/dh and dU/dh (layer thickness), not
%    density. Those are now returned as dcdh/dudh, and dcdrho/dudrho are
%    true density kernels computed with calc_density_kernel96.m.
%  - Spherical-earth correction of the phase-velocity kernels: the kernels
%    are computed by srfdrr96/srfdrl96 for the earth-flattened model and
%    must be multiplied by dv_flat/dv_sph to apply to the spherical model.
%    Because only group velocity is requested here, srfdrr96/srfdrl96
%    apply that mapping to the group-velocity kernels only and output the
%    phase-velocity kernels dc/dVs and dc/dVp unmapped. They are now
%    multiplied by the flattening factor 2a/(r_top+r_bot) (a = 6370 km),
%    which can reach several percent in the upper mantle (2.4% at 150 km).
%    This was verified against finite differences and against the kernels
%    of sregn96/slegn96 -DER.
%

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

if size(model,2)>4 % if Qp and Qs proivded, calculate the Q-corrected kernels

    [kernel, dispersion] = calc_anelastic_kernel96(model, vec_T, wavetype, ifnorm, ifplot, nmode, fref);
    dcdvs=kernel.dcdvs;
    dcdvp=kernel.dcdvp;
    dudvs=kernel.dudvs;
    dudvp=kernel.dudvp;
    dcdrho=kernel.dcdrho;
    dudrho=kernel.dudrho;
    dcdh=kernel.dcdh;
    dudh=kernel.dudh;
    z_mid=kernel.z;
    phv=dispersion.phv;
    grv=dispersion.grv;

else % Otherwise, calculate non-Q-corrected kernels

    % timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
    timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

    [nlayer,temp]=size(model);
    Nper = length(vec_T);
    dcdvs=nan*ones(nlayer,Nper);
    dcdvp=nan*ones(nlayer,Nper);
    dudvs=nan*ones(nlayer,Nper);
    dudvp=nan*ones(nlayer,Nper);
    dcdh=nan*ones(nlayer,Nper);
    dudh=nan*ones(nlayer,Nper);

    % Earth-flattening factor dv_flat/dv_sph for each layer (subroutine
    % sphere of srfdis96: a = 6370 km, mid-layer rule, 1 km for the
    % half-space). writemod_surf96 writes SPHERICAL EARTH models.
    a_earth = 6370;
    ztop = [0; cumsum(model(1:end-1,1))];
    hh = model(:,1); hh(end) = 1;
    btp = 2*a_earth ./ ((a_earth - ztop) + (a_earth - ztop - hh));

    % Density kernels and elastic phase/group velocity (sregn96/slegn96)
    if(wavetype=='J')
        dcdrho=nan*ones(nlayer,Nper);
        dudrho=nan*ones(nlayer,Nper);
        phv=nan*ones(Nper,1);
        grv=nan*ones(Nper,1);
    else
        [dcdrho, dudrho, phv, grv] = calc_density_kernel96(model, vec_T, wavetype, nmode);
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
            dcdh = [];
            dudh = [];
            z_mid = [];
            return
        end
        % phase-velocity kernels: apply the spherical-earth mapping (see header)
        dcdvs(:,ip)=kernel(:,1).*btp;
        dcdvp(:,ip)=kernel(:,2).*btp;
        dudvs(:,ip)=kernel(:,3);
        dudvp(:,ip)=kernel(:,4);
        dcdh(:,ip)=kernel(:,5);
        dudh(:,ip)=kernel(:,6);

        if(ifnorm)
            dcdvs(:,ip) = dcdvs(:,ip)./model(:,1);
            dcdvp(:,ip) = dcdvp(:,ip)./model(:,1);
            dudvs(:,ip) = dudvs(:,ip)./model(:,1);
            dudvp(:,ip) = dudvp(:,ip)./model(:,1);
            dcdrho(:,ip) = dcdrho(:,ip)./model(:,1);
            dudrho(:,ip) = dudrho(:,ip)./model(:,1);
            dcdh(:,ip) = dcdh(:,ip)./model(:,1);
            dudh(:,ip) = dudh(:,ip)./model(:,1);
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
        dcdh(isinf(dcdh(:,ip)),ip) = 0;
        dudh(isinf(dudh(:,ip)),ip) = 0;

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
end

return
end
