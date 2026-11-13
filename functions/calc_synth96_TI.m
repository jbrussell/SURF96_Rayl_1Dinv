% Synthetic seismograms by surface-wave mode summation for a transversely
% isotropic (TI) model, using the CPS programs tprep96, tdisp96,
% tlegn96/tregn96, tpulse96, fmech96 and f96tosac
%
%   syn = calc_synth96_TI(model, dist, par, src)
%
% model: TI model [H Vpv Vsv rho Vph Vsh eta] or [... Qp Qs]
%        (see writemod_surf96_TI.m; Qp = 1000, Qs = 200 if not given)
% dist:  epicentral distances (km) [vector]
% par:   structure of options (defaults in brackets)
%   .hs        source depth (km, below the sea surface) [required]
%   .hr        receiver depth (km) [0]. A receiver exactly at the bottom of
%              a water layer (a sea-floor station) gives zero Love-wave
%              motion in CPS, so it is moved 0.01 km into the solid.
%   .dt        sampling interval (s) [1]; Nyquist = 1/(2 dt)
%   .npts      number of samples [4096]; the synthetics are periodic with
%              period npts*dt: make it longer than the whole signal
%   .t0, .vred first sample at t0 + dist/vred (vred = 0: at t0) [0, 0]
%   .fmax      highest frequency computed (Hz) [Nyquist]; lower values
%              save time (TI Rayleigh dispersion is slow)
%   .nmod      maximum number of modes computed [20]
%   .wavetype  'R', 'L' or 'RL' [ 'RL' ]
%   .modesel   cell array of mode selections to synthesize, each 'all'
%              (sum of all computed modes) or a mode number (0 = fund.,
%              1 = first overtone, ...) [{'all'}]
%   .pulse     source time function: 't' triangle (base 2*L*dt),
%              'p' parabolic (base 4*L*dt), 'i' impulse ['t']
%   .L         pulse duration parameter (samples) [1]
%   .output    'D' displacement (cm), 'V' velocity (cm/s), 'A' acceleration
%              (cm/s^2) ['D']
%   .fref      reference frequency (Hz) of the Q model [1]
% src:   source structure
%   .az        source-to-station azimuth (deg, from north) [required]
%   and one of
%   .strike, .dip, .rake (deg) with .Mw or .M0 (dyne-cm)   double couple
%   .mt = [Mxx Myy Mzz Mxy Mxz Myz] (dyne-cm, x = N, y = E, z = down)
%   .force = [Fn Fe Fdown] (dyne)                           point force
%
% syn: structure
%   .t          time (s) [npts x ndist]
%   .Z, .R, .T  vertical (up), radial and transverse ground motion
%               [npts x ndist x nsel] (in the units of par.output)
%   .modesel    the mode selections (as par.modesel)
%   .dist, .az, .par
%
% Mode numbers are those of tdisp96 (ordered by phase velocity). Unlike
% calc_der96_TI.m, tpulse96 uses the tdisp96 roots without rejecting the
% occasional spurious or duplicate root (see calc_der96_TI.m), which can
% shift the numbering of the modes above it at some frequencies.
% Body waves (modes faster than the half-space S velocity) are not
% included.
%
% The CPS programs tpulse96 (VOLIII), fmech96 and f96tosac (VOLV) are not in
% bin_v3.30/; build them from src.zip (PROGRAMS.330) like the other CPS
% programs and put them on the PATH (or in the CPS_BIN directory).
%
% jbrussell 11/2026
%
function syn = calc_synth96_TI(model, dist, par, src)

% defaults
def = struct('hr',0,'dt',1,'npts',4096,'t0',0,'vred',0,'fmax',[],'nmod',20, ...
    'wavetype','RL','modesel',{{'all'}},'pulse','t','L',1,'output','D','fref',1);
fn = fieldnames(def);
for k = 1:length(fn)
    if ~isfield(par,fn{k}) || isempty(par.(fn{k}))
        par.(fn{k}) = def.(fn{k});
    end
end
if ~isfield(par,'hs'); error('calc_synth96_TI: par.hs (source depth) is required'); end
if ~isfield(src,'az'); error('calc_synth96_TI: src.az (source-to-station azimuth) is required'); end
if ~iscell(par.modesel); par.modesel = {par.modesel}; end
dist = dist(:)';
wt = upper(par.wavetype);

% CPS programs: full paths (the PATH may hold relative directories such as
% ./bin_v3.30/, and the programs are run in a temporary directory)
progs = {'tprep96','tdisp96','tlegn96','tregn96','tpulse96','fmech96','f96tosac'};
exe = struct();
for k = 1:length(progs)
    [st, p] = system(['command -v ',progs{k}]);
    p = strtrim(p);
    if st ~= 0 || isempty(p)
        error(['calc_synth96_TI: ',progs{k},' is not on the PATH. Build tpulse96 (VOLIII), ', ...
            'fmech96 and f96tosac (VOLV) from src.zip (PROGRAMS.330) and add them to the PATH.']);
    end
    if p(1) ~= '/'
        p = fullfile(pwd, p);
    end
    exe.(progs{k}) = ['"',p,'"'];
end

% receiver exactly at the bottom of a water layer -> just inside the solid
ztop = [0; cumsum(model(1:end-1,1))];
zbot = ztop + model(:,1);
ifl = find(model(:,3) <= 0);
hr = par.hr;
if any(abs(zbot(ifl) - hr) < 1e-6)
    hr = hr + 0.01;
end

% source mechanism arguments for fmech96
if isfield(src,'strike')
    mech = sprintf(' -S %.6g -D %.6g -R %.6g', src.strike, src.dip, src.rake);
    if isfield(src,'Mw')
        mech = [mech, sprintf(' -MW %.6g', src.Mw)];
    elseif isfield(src,'M0')
        mech = [mech, sprintf(' -M0 %.6g', src.M0)];
    end
    gfopt = ' -EQEX';
elseif isfield(src,'mt')
    m = src.mt;
    mech = sprintf(' -XX %.6g -YY %.6g -ZZ %.6g -XY %.6g -XZ %.6g -YZ %.6g', m);
    gfopt = ' -EQEX';
elseif isfield(src,'force')
    f = src.force;
    mech = sprintf(' -fx %.6g -fy %.6g -fz %.6g', f);
    gfopt = ' -EXF';
else
    error('calc_synth96_TI: src needs strike/dip/rake, mt or force');
end
mech = [mech, sprintf(' -A %.6g -ROT', src.az)];

switch lower(par.pulse)
    case 't'; pulse = sprintf(' -t -l %d', par.L);
    case 'p'; pulse = sprintf(' -p -l %d', par.L);
    case 'i'; pulse = ' -i';
    otherwise; error('calc_synth96_TI: par.pulse must be ''t'', ''p'' or ''i''');
end
outopt = [' -', upper(par.output)];

% work in a temporary directory (the CPS programs are run there; MATLAB's
% current directory is not changed, so relative paths stay valid)
tmpd = tempname;
mkdir(tmpd);
cleanupObj = onCleanup(@() cleanup_tmp(tmpd));
runcps = @(cmd) system(['cd "', tmpd, '" && ', cmd]);
writemod_surf96_TI(model, fullfile(tmpd,'model.mod'), par.fref);

fid = fopen(fullfile(tmpd,'dfile'),'w');
for k = 1:length(dist)
    fprintf(fid, '%.4f %.6g %d %.6g %.6g\n', dist(k), par.dt, par.npts, par.t0, par.vred);
end
fclose(fid);

% dispersion and eigenfunctions
prepopt = sprintf(' -HS %.6g -HR %.6g -NMOD %d', par.hs, hr, par.nmod);
if any(wt == 'L'); prepopt = [prepopt, ' -L']; end
if any(wt == 'R'); prepopt = [prepopt, ' -R']; end
if ~isempty(par.fmax); prepopt = [prepopt, sprintf(' -FMAX %.6g', par.fmax)]; end
[st, log] = runcps([exe.tprep96,' -M model.mod -d dfile', prepopt]);
if st ~= 0; error('calc_synth96_TI: tprep96 failed:\n%s', log); end
% tdisp96 limits the phase-velocity search to the largest SV velocity;
% raise it to the largest SH velocity for Love waves (as in calc_der96_TI)
tdispopt = '';
if any(wt == 'L') && any(model(:,6) > model(:,3))
    a_e = 6371;
    hh = model(:,1); hh(end) = 1;
    fac = 2*a_e ./ ((a_e - ztop) + (a_e - ztop - hh));
    tdispopt = sprintf(' -cmax %.5f', max(model(:,6).*fac));
end
[st, log] = runcps([exe.tdisp96, tdispopt]);
if st ~= 0; error('calc_synth96_TI: tdisp96 failed:\n%s', log); end
if any(wt == 'L'); [st, log] = runcps(exe.tlegn96); if st ~= 0; error('calc_synth96_TI: tlegn96 failed:\n%s', log); end; end
if any(wt == 'R'); [st, log] = runcps(exe.tregn96); if st ~= 0; error('calc_synth96_TI: tregn96 failed:\n%s', log); end; end

% mode summation, source mechanism and conversion to SAC
nd = length(dist); nsel = length(par.modesel);
syn.Z = zeros(par.npts, nd, nsel); syn.R = syn.Z; syn.T = syn.Z;
syn.t = zeros(par.npts, nd);
comp = {'Z','R','T'};
for is = 1:nsel
    ms = par.modesel{is};
    if ischar(ms)
        if ~strcmpi(ms,'all'); error('calc_synth96_TI: modesel entries must be ''all'' or a mode number'); end
        mopt = '';
    else
        mopt = sprintf(' -M %d', ms);
    end
    delete_if(fullfile(tmpd,'*.sac'));
    [st, log] = runcps([exe.tpulse96,' -d dfile', pulse, outopt, gfopt, mopt, ' > gf.f96']);
    if st ~= 0; error('calc_synth96_TI: tpulse96 failed:\n%s', log); end
    [st, log] = runcps([exe.fmech96, mech, ' < gf.f96 > syn.f96']);
    if st ~= 0; error('calc_synth96_TI: fmech96 failed:\n%s', log); end
    [st, log] = runcps([exe.f96tosac,' -B syn.f96']);
    if st ~= 0; error('calc_synth96_TI: f96tosac failed:\n%s', log); end
    for k = 1:nd
        for ic = 1:3
            f = fullfile(tmpd, sprintf('B%03d%02d%s00.sac', k, ic, comp{ic}));
            if ~exist(f,'file'); continue; end
            [x, hdr] = readsac(f);
            n = min(length(x), par.npts);
            syn.(comp{ic})(1:n,k,is) = x(1:n);
            syn.t(:,k) = hdr.b + (0:par.npts-1)'*hdr.delta;
        end
    end
end
syn.modesel = par.modesel;
syn.dist = dist;
syn.az = src.az;
par.hr_used = hr;
syn.par = par;

end

function [x, hdr] = readsac(f)
% read a little-endian binary SAC file
fid = fopen(f,'r','ieee-le');
h = fread(fid,70,'float32');
hi = fread(fid,40,'int32');
fseek(fid,632,'bof');
x = fread(fid,hi(10),'float32');
fclose(fid);
hdr.delta = h(1); hdr.b = h(6); hdr.npts = hi(10);
end

function delete_if(pattern)
d = dir(pattern);
for k = 1:length(d); delete(fullfile(fileparts(pattern), d(k).name)); end
end

function cleanup_tmp(tmpd)
if exist(tmpd,'dir')
    rmdir(tmpd,'s');
end
end
