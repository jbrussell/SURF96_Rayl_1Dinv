function [dcdrho, dudrho, phv, grv] = calc_density_kernel96(model, vec_T, wavetype, varargin)
% Calculate phase- and group-velocity density (rho) sensitivity kernels for
% an isotropic model. srfker96 does not output density kernels, so they are
% computed here with sprep96 / sdisp96 / sregn96 (slegn96) -DER.
%
% [dcdrho, dudrho, phv, grv] = calc_density_kernel96(model, vec_T, wavetype);
% [dcdrho, dudrho, phv, grv] = calc_density_kernel96(model, vec_T, wavetype, nmode);
%
% model:    [thickness, vp, vs, rho] (any Qp, Qs columns are ignored: elastic)
% vec_T:    periods (s)
% wavetype: 'R' or 'L'
% nmode:    0 = fundamental mode, 1 = 1st overtone, ... (OPTIONAL, default 0)
%
% dcdrho, dudrho: [nlayer x nperiod] layer kernels dc/drho and dU/drho
%                 (velocities held fixed), not normalized by thickness
% phv, grv:       elastic phase and group velocity at vec_T (column vectors)
% Periods at which the mode does not exist are returned as NaN.
%
% dc/drho is read from the binary s?egn96.der file (read_der96_TI.m reads
% both the isotropic and TI files) at full single precision. dU/drho is
% computed from dc/drho following Rodi et al. (1975):
%     dU/dm = (U/c)(2 - U/c) dc/dm + (U^2 w / c^2) d(dc/dm)/dw
% with a centered finite difference over w*(1 +/- dw_frac), dw_frac = 2e-2.
%
% jbrussell 11/2026
%

if nargin >= 4
    nmode = varargin{1};
else
    nmode = 0;
end
wavetype = upper(wavetype);

% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

dw_frac = 2e-2;
vec_T = vec_T(:)';
Nper = length(vec_T);
nlayer = size(model,1);
T_all = [vec_T, vec_T./(1+dw_frac), vec_T./(1-dw_frac)];

files = {'periods.txt','start.mod','sdisp96.dat','sdisp96.ray','sdisp96.lov', ...
    'sregn96.der','slegn96.der'};
cleanup(files);

% elastic model (Q^-1 = 0)
writemod_surf96([model(:,1:4), zeros(nlayer,2)],'start.mod');
fid = fopen('periods.txt','w');
fprintf(fid,'%.8f\n',T_all);
fclose(fid);
[~, log] = system([timeoutstr,'sprep96 -M start.mod -NMOD ',num2str(nmode+1),' -',wavetype,' -PARR periods.txt']);
[~, log] = system([timeoutstr,'sdisp96']);
[~, log] = system([timeoutstr,'s',lower(wavetype),'egn96 -DER']);

dcdrho = nan(nlayer,Nper);
dudrho = nan(nlayer,Nper);
phv = nan(Nper,1);
grv = nan(Nper,1);
derfile = ['s',lower(wavetype),'egn96.der'];
if ~exist(derfile,'file')
    cleanup(files);
    warning('calc_density_kernel96: %s was not created', derfile);
    return
end
eig = read_der96_TI(derfile, wavetype, nmode);
cleanup(files);

for ip = 1:Nper
    i0 = findper(eig.periods, vec_T(ip));
    ipl = findper(eig.periods, vec_T(ip)/(1+dw_frac)); % w(1+dw_frac)
    imn = findper(eig.periods, vec_T(ip)/(1-dw_frac)); % w(1-dw_frac)
    if i0 == 0 || isnan(eig.phv(i0))
        continue
    end
    c = eig.phv(i0); U = eig.grv(i0);
    phv(ip) = c; grv(ip) = U;
    K = eig.dcdrho(:,i0);
    dcdrho(:,ip) = K;
    okp = ipl > 0 && ~isnan(eig.phv(ipl));
    okm = imn > 0 && ~isnan(eig.phv(imn));
    if okp && okm
        wdKdw = (eig.dcdrho(:,ipl) - eig.dcdrho(:,imn)) / (2*dw_frac);
    elseif okp
        wdKdw = (eig.dcdrho(:,ipl) - K) / dw_frac;
    elseif okm
        wdKdw = (K - eig.dcdrho(:,imn)) / dw_frac;
    else
        continue
    end
    dudrho(:,ip) = (U/c)*(2-U/c)*K + (U^2/c^2)*wdKdw;
end

end

function I = findper(periods, T)
[dT, I] = min(abs(periods - T));
if isempty(I) || dT > 1e-5*T
    I = 0;
end
end

function cleanup(files)
for ii = 1:length(files)
    if exist(files{ii},'file'); delete(files{ii}); end
end
end
