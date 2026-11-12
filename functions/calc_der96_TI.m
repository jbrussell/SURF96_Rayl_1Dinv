% Run the CPS transversely isotropic (TI) surface-wave programs and return
% validated eigenfunctions and phase-velocity kernels for one mode branch.
% This is the engine behind calc_eigenfunctions96_TI.m, calc_kernel96_TI.m,
% calc_amplification_ellipticity_gamma_TI.m and dispR/dispL_surf96_TI.m.
%
% [eig, eig_p, eig_m, info] = calc_der96_TI(vec_T, model, wavetype, nmode, fref)
% [eig, eig_p, eig_m, info] = calc_der96_TI(vec_T, model, wavetype, nmode, fref, opts)
%
% 	vec_T:    periods (s)
% 	model:    TI model [H Vpv Vsv rho Vph Vsh eta] or [... Qp Qs]
%             (see writemod_surf96_TI.m; Qp=1000, Qs=200 if not given)
%   wavetype: 'R' or 'L'
%   nmode:    mode branch (0 = fund.), or 'maxamp' (see below)
%   fref:     reference frequency of Q model (Hz)
%   opts:     (OPTIONAL) structure with any of the fields
%     .nmodes  number of modes (roots) tdisp96 computes at each period
%              (default nmode+1, i.e. 1 for the fundamental mode; N+1 for
%              nmode = 'maxampN'). Periods
%              at which no valid root of branch nmode is found are computed
%              again once with nmodes+2 (spurious roots can occupy a slot).
%     .nfreq   frequencies computed per period: 3 = w, w(1+dw_frac),
%              w(1-dw_frac) (default; needed for the group-velocity kernels
%              of calc_kernel96_TI and to rebuild missed roots), or 2 = w and
%              w(1+dw_frac) (enough for phase-velocity kernels and root
%              validation; eig_m is then NaN and missed roots are not rebuilt)
%     .c_guess expected phase velocity at each period [1 x Nper] (km/s),
%              e.g. a linearized prediction in an inversion. Only w is
%              computed and the root closest to c_guess (within 2%) is
%              selected; periods without such a root are computed again
%              with nfreq = 2 and the full root validation. eig_p and eig_m
%              are NaN. (With nmode = 'maxamp' this follows the branch of
%              the guess; eig.mode is then the position of the selected
%              root among the roots found.)
%
% eig:   output of read_der96_TI.m at the requested periods (in the order
%        requested; NaN where the mode does not exist)
% eig_p, eig_m: the same mode branch at angular frequency w*(1+dw_frac) and
%        w*(1-dw_frac) (used for group-velocity kernels, calc_kernel96_TI.m)
% info:  .dw_frac, .U_tol
%        .nspurious  number of roots discarded at each requested period
%                    (spurious or duplicate roots below the selected mode)
%        .repaired   true where the selected root was missing at the
%                    requested period and was rebuilt from the neighbors
%        .raw        all roots found by tdisp96 at the requested periods
%                    (.periods, .phv, .grv)
%        .nretry     number of periods recomputed with more modes
%        .nfallback  number of periods where no root matched c_guess
%
% ROOT VALIDATION
% tdisp96 searches for sign changes of the dispersion function in steps of
% wavenumber and refines each with 15 bisections. For TI models with a
% water layer over slow sediments it occasionally (in tests ~1-5% of
% periods; the isotropic sdisp96 and surf96 did not do this for the same
% models)
%   - reports roots that are not modes of the model (nearly non-dispersive,
%     typically near the slowest S velocity),
%   - reports the same root twice, or
%   - misses a root.
% Any of these changes the numbering of the modes above it. The modes are
% therefore computed at three frequencies, w*(1-dw_frac), w and
% w*(1+dw_frac), and followed as branches: two roots at adjacent
% frequencies belong to the same branch if the secant group velocity
% dw/dk between them matches the mean of their group velocities from the
% energy integrals (tregn96/tlegn96) to within U_tol. For genuine roots the
% misfit is ~1e-6 to 3e-3 (limited by the precision of the tdisp96 roots);
% spurious roots are off by several percent. A branch is accepted if it is
% found consistently at at least two of the three frequencies; isolated
% roots and duplicates are discarded. If an accepted branch was missed at
% the requested frequency, its values there are rebuilt as the average of
% the two neighbors (second-order accurate; eigenfunctions are first
% brought to a common amplitude). The mode branch nmode is then selected
% from the ordered list of accepted branches. Roots at the (flattened)
% half-space S velocity, the upper end of the tdisp96 search range, where
% it can report non-dispersive roots, are discarded.
% A TI:spuriousRoot warning is issued whenever roots are discarded or
% repaired; disable it with warning('off','TI:spuriousRoot').
%
% MAXIMUM SURFACE AMPLITUDE BRANCH (nmode = 'maxamp' or 'maxampN')
% Modes are numbered by phase velocity. In models with a low-velocity
% waveguide (e.g. a thin oceanic crust over a fast lid and an LVZ, such as
% PA5) the slowest mode at short periods can be trapped at depth and
% barely reach the surface, while the wave that dominates the surface
% records is mode 1 (mode osculation). With nmode = 'maxamp', at each
% period the branch among modes 0..N (N = 2 by default, 'maxampN' for
% another N) with the largest surface amplitude for a source and a
% receiver at the top of the solid (sea floor) is selected:
%   Love:     AL * ut(top)^2      Rayleigh: AR * uz(top)^2
% (AL, AR = 1/(2 c U I0) are the CPS energy factors; the product does not
% depend on the arbitrary eigenfunction amplitude). The neighbors at
% w(1+-dw_frac) are taken from the same branch. eig.mode then holds the
% selected mode number at each period. The selection can change from one
% model to the next (e.g. between inversion iterations), so the phase
% velocity is not a smooth function of the model where two branches have
% similar amplitudes.
%
% Q
% CPS is always run with Q^-1 = 0. tregn96/tlegn96 correct the phase
% velocity for causal (physical) dispersion but not the group velocity,
% which would bias the root validation. The attenuation coefficient and the
% causal phase velocity are computed here afterwards with the formulas of
% tregn96/tlegn96 (subroutine gammap) from the elastic kernels; they agree
% with the CPS values to ~1e-7. eig.phv_elastic holds the elastic phase
% velocity. As in CPS, group velocity, eigenfunctions and kernels are the
% elastic values.
%
% COST
% Nearly all of the time is spent in tdisp96 (Rayleigh waves especially).
% The cost is proportional to nfreq x number of periods, and grows with
% nmodes (about 1.6x for 3 modes instead of 1).
%
% ONE PERIOD AT A TIME
% When tdisp96 is given several periods, the number of roots it accepts at
% each period is capped by the number it found at the previous (shorter)
% period. With a spurious root this drops genuine modes at all longer
% periods, so tprep96/tdisp96/tregn96 are run separately for every period
% (nfreq per requested period), each with the 20 s gtimeout used by the
% isotropic wrappers. tdisp96 is considerably slower than sdisp96 for
% Rayleigh waves (~0.1-0.2 s per period for a ~40 layer model).
%
% jbrussell 11/2026
%
function [eig, eig_p, eig_m, info] = calc_der96_TI(vec_T,model,wavetype,nmode,fref,varargin)

dw_frac = 2e-2;     % fractional frequency step for neighbors
U_tol = 1.2e-2;     % relative group-velocity consistency tolerance
dup_tol = 1e-5;     % relative phase-velocity difference below which two roots are duplicates
NRETRY = 2;         % extra modes computed when a period has to be recomputed
c_tol = 2e-2;       % max relative difference from c_guess

opts = struct();
if ~isempty(varargin) && isstruct(varargin{1})
    opts = varargin{1};
end
% mode selection: fixed branch number, or the branch with the largest
% surface amplitude among modes 0..nmax_amp
maxamp = ischar(nmode);
nmax_amp = NaN;
if maxamp
    if ~strncmpi(nmode,'maxamp',6)
        error('calc_der96_TI: nmode must be a mode number or ''maxamp''');
    end
    nmax_amp = 2;
    if length(nmode) > 6
        nmax_amp = str2double(nmode(7:end));
    end
    nsearch = nmax_amp; % highest mode number needed
else
    nsearch = nmode;
end
isol = find(model(:,3) > 0, 1); % top of the solid (sea floor) for 'maxamp'

nmodes = nsearch + 1;
if isfield(opts,'nmodes') && ~isempty(opts.nmodes)
    nmodes = opts.nmodes;
end
if nmodes < nsearch + 1
    error('calc_der96_TI: opts.nmodes (%d) must be at least %d for this nmode', nmodes, nsearch+1);
end
nfreq = 3;
if isfield(opts,'nfreq') && ~isempty(opts.nfreq)
    nfreq = opts.nfreq;
end
if ~any(nfreq == [2 3])
    error('calc_der96_TI: opts.nfreq must be 2 or 3');
end
c_guess = [];
if isfield(opts,'c_guess') && ~isempty(opts.c_guess)
    c_guess = opts.c_guess(:)';
end

wavetype = upper(wavetype);
vec_T = vec_T(:)';
Nper = length(vec_T);
if ~isempty(c_guess) && length(c_guess) ~= Nper
    error('calc_der96_TI: opts.c_guess must have one value per period');
end

cleanup_cps96_TI(); % clean up

% Warn about a known CPS issue (see check_halfspace_eta_TI.m)
check_halfspace_eta_TI(model,wavetype);

% Q model (as in writemod_surf96_TI: Qp=1000, Qs=200 if not given) and
% CPS convention for Q^-1 (Q>1 is converted to 1/Q; Q=0 means elastic)
nlay = size(model,1);
if size(model,2) >= 9
    qa = model(:,8); qb = model(:,9);
else
    qa = 1000*ones(nlay,1); qb = 200*ones(nlay,1);
end
qa(qa>1) = 1./qa(qa>1);
qb(qb>1) = 1./qb(qb>1);
qb(model(:,3)<=0.01) = 0; % fluid

% Make model file. CPS is always run with Q^-1 = 0: tregn96/tlegn96 apply
% the causal (physical dispersion) correction to the phase velocity but
% not to the group velocity, which would break the root validation below.
% The Q correction is applied afterwards (see APPLY Q below).
writemod_surf96_TI([model(:,1:7), zeros(nlay,2)],'start.mod',fref);

% tdisp96 limits the phase-velocity search to the largest SV velocity (TL)
% of the flattened model, also for Love waves. When Vsh > Vsv this misses
% Love modes faster than that, so for Love waves the upper limit is set to
% the largest flattened SH velocity instead (-cmax).
tdispopt = '';
if strcmp(wavetype,'L') && any(model(:,6) > model(:,3))
    a_e = 6371;
    ztop = [0; cumsum(model(1:end-1,1))];
    hh = model(:,1); hh(end) = 1; % tdisp96 uses 1 km for the half-space
    fac = 2*a_e ./ ((a_e - ztop) + (a_e - ztop - hh));
    tdispopt = sprintf(' -cmax %.5f', max(model(:,6).*fac));
end

% half-space shear velocity (SH for Love, SV for Rayleigh), mapped to the
% flattened model used by tdisp96 (writemod_surf96_TI writes SPHERICAL
% EARTH models; tdisp96 uses radius 6371 km and the layer mid-point rule).
% tdisp96 can report non-dispersive roots right at this velocity (the upper
% end of its search range)
if strcmp(wavetype,'L')
    vs_hs = model(end,6);
else
    vs_hs = model(end,3);
end
a_earth = 6371;
z_hs = sum(model(1:end-1,1));
vs_hs = vs_hs * 2*a_earth/(2*(a_earth - z_hs) - 1);

% Frequency factors: 1 = w, 2 = w(1+dw_frac), 3 = w(1-dw_frac)
ffac = [1, 1+dw_frac, 1-dw_frac];
R = cell(Nper,3);       % roots at each period and frequency (see run_one)
nmod_used = zeros(1,Nper);
ntimeout = 0;

% ---- First pass ----
if ~isempty(c_guess)
    nf_first = 1;
else
    nf_first = nfreq;
end
for ii = 1:Nper
    for f = 1:nf_first
        [R{ii,f}, nto] = run_one(vec_T(ii)/ffac(f), nmodes, wavetype, tdispopt);
        ntimeout = ntimeout + nto;
    end
    nmod_used(ii) = nmodes;
end

% ---- Root selection ----
sel = repmat(struct('jz',NaN,'jp',NaN,'jm',NaN,'repaired',false,'nspur',0,'mode',NaN), 1, Nper);
nfallback = 0;
nretry = 0;
for ii = 1:Nper
    ok = false;
    if ~isempty(c_guess)
        % pick the root closest to the expected phase velocity
        [cz, ~, jz, ndz] = getroots(R{ii,1}, dup_tol, vs_hs);
        if ~isempty(cz) && ~isnan(c_guess(ii))
            [dmin, k] = min(abs(cz - c_guess(ii)));
            if dmin <= c_tol*c_guess(ii)
                sel(ii).jz = jz(k); sel(ii).nspur = ndz;
                if maxamp
                    sel(ii).mode = k - 1;
                else
                    sel(ii).mode = nmode;
                end
                ok = true;
            end
        end
        if ~ok
            % fall back to the full validation with two frequencies
            nfallback = nfallback + 1;
            [R{ii,2}, nto] = run_one(vec_T(ii)/ffac(2), nmodes, wavetype, tdispopt);
            ntimeout = ntimeout + nto;
            sel(ii) = select_branch(R(ii,:), 2, vec_T(ii), nmode, dw_frac, U_tol, dup_tol, vs_hs, maxamp, nmax_amp, wavetype, isol);
            ok = ~isnan(sel(ii).jz) || sel(ii).repaired;
        end
        nf_ii = 2;
    else
        sel(ii) = select_branch(R(ii,:), nfreq, vec_T(ii), nmode, dw_frac, U_tol, dup_tol, vs_hs, maxamp, nmax_amp, wavetype, isol);
        ok = ~isnan(sel(ii).jz) || sel(ii).repaired;
        nf_ii = nfreq;
    end
    if ~ok && nmod_used(ii) < nmodes + NRETRY
        % recompute this period with more modes (a spurious root may have
        % taken the place of the requested one)
        nretry = nretry + 1;
        for f = 1:nf_ii
            [R{ii,f}, nto] = run_one(vec_T(ii)/ffac(f), nmodes + NRETRY, wavetype, tdispopt);
            ntimeout = ntimeout + nto;
        end
        nmod_used(ii) = nmodes + NRETRY;
        sel(ii) = select_branch(R(ii,:), nf_ii, vec_T(ii), nmode, dw_frac, U_tol, dup_tol, vs_hs, maxamp, nmax_amp, wavetype, isol);
    end
end
if ntimeout > 0
    warning('TI:timeout', 'calc_der96_TI: a CPS program was stopped by gtimeout %d time(s); results there are missing', ntimeout);
end

% ---- Assemble the output ----
tmpl = [];
for ii = 1:Nper
    for f = 1:3
        if ~isempty(R{ii,f}) && ~isempty(R{ii,f}.eig); tmpl = R{ii,f}.eig; break; end
    end
    if ~isempty(tmpl); break; end
end
eig   = pickroots(R, 1, [sel.jz], tmpl, Nper);
eig_p = pickroots(R, 2, [sel.jp], tmpl, Nper);
eig_m = pickroots(R, 3, [sel.jm], tmpl, Nper);
eig.periods = vec_T;
eig_p.periods = vec_T/(1+dw_frac);
eig_m.periods = vec_T/(1-dw_frac);
if maxamp
    modesel = [sel.mode];
    eig.mode = modesel; eig_p.mode = modesel; eig_m.mode = modesel;
else
    eig.mode = nmode; eig_p.mode = nmode; eig_m.mode = nmode;
end
repaired = [sel.repaired];
nspurious = [sel.nspur];

% Rebuild roots missed at the requested frequency from the two neighbors
if any(repaired)
    if strcmp(wavetype,'L')
        efun = {'ut','tt'}; eref = 'ut';
    else
        efun = {'ur','tr','uz','tz'}; eref = 'uz';
    end
    flds = fieldnames(eig);
    for ii = find(repaired)
        % bring the eigenfunction at w(1-dw_frac) to the amplitude of the
        % one at w(1+dw_frac) (amplitude is arbitrary; A scales as 1/amp^2)
        a = eig_p.(eref)(:,ii); bm = eig_m.(eref)(:,ii);
        s = (a'*bm)/(bm'*bm);
        em = struct();
        for jj = 1:length(flds)
            if any(strcmp(flds{jj},{'periods','mode'})); continue; end
            em.(flds{jj}) = eig_m.(flds{jj})(:,ii);
        end
        for jj = 1:length(efun)
            em.(efun{jj}) = em.(efun{jj})*s;
        end
        em.amp = em.amp/s^2;
        for jj = 1:length(flds)
            if any(strcmp(flds{jj},{'periods','mode'})); continue; end
            eig.(flds{jj})(:,ii) = (eig_p.(flds{jj})(:,ii) + em.(flds{jj}))/2;
        end
    end
end

% APPLY Q: attenuation coefficient and causal phase velocity, following
% subroutine gammap of tregn96/tlegn96:
%   gamma = w/(2 c0^2) sum_i [ (dc/dVsv Vsv + dc/dVsh Vsh) Qs^-1
%                            + (dc/dVpv Vpv + dc/dVph Vph) Qp^-1 ]
%   c = c0 + sum_i [ ... ] ln(f/fref)/pi
% (eigenfunctions, group velocity and kernels are the elastic values, as in
% CPS)
eig   = applyQ(eig, model, qa, qb, fref);
eig_p = applyQ(eig_p, model, qa, qb, fref);
eig_m = applyQ(eig_m, model, qa, qb, fref);

% table of the roots found at the requested periods
maxr = 0;
for ii = 1:Nper
    if ~isempty(R{ii,1}); maxr = max(maxr, length(R{ii,1}.c)); end
end
raw.periods = vec_T;
raw.phv = nan(Nper, maxr); raw.grv = nan(Nper, maxr);
for ii = 1:Nper
    if ~isempty(R{ii,1})
        n = length(R{ii,1}.c);
        raw.phv(ii,1:n) = R{ii,1}.c; raw.grv(ii,1:n) = R{ii,1}.u;
    end
end

info.dw_frac = dw_frac;
info.U_tol = U_tol;
info.nspurious = nspurious;
info.repaired = repaired;
info.nretry = nretry;
info.nfallback = nfallback;
info.raw = raw;
if any(nspurious > 0) || any(repaired)
    warning('TI:spuriousRoot', ['calc_der96_TI: tdisp96 root problems at period(s) %s s: discarded %d spurious/duplicate root(s), ', ...
        'rebuilt %d missed root(s) from neighboring frequencies (see calc_der96_TI.m). ', ...
        'Disable this message with warning(''off'',''TI:spuriousRoot'').'], ...
        mat2str(vec_T(nspurious>0 | repaired),4), sum(nspurious), sum(repaired));
end

cleanup_cps96_TI();

end

function [Rk, ntimeout] = run_one(T, nmod, wavetype, tdispopt)
% Run tprep96, tdisp96 and tregn96/tlegn96 -DER for one period (the model
% file start.mod must exist). Rk.c, Rk.u: phase and group velocity of all
% roots; Rk.eig: read_der96_TI output with one column per root (empty if
% none).
% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time
derfile = ['t',lower(wavetype),'egn96.der'];
if exist(derfile,'file'); delete(derfile); end
fid = fopen('periods.txt','w');
fprintf(fid,'%.8f\n',T);
fclose(fid);
% tprep96: dispersion control file (tdisp96.dat)
[st1, log] = system([timeoutstr,'tprep96 -M start.mod -NMOD ',num2str(nmod),' -',wavetype,' -PARR periods.txt']);
% tdisp96: dispersion binary file (tdisp96.ray or tdisp96.lov)
[st2, log] = system([timeoutstr,'tdisp96',tdispopt]);
% tregn96/tlegn96 -DER: eigenfunctions and partial derivatives
[st3, log] = system([timeoutstr,'t',lower(wavetype),'egn96 -DER']);
ntimeout = any([st1 st2 st3] == 124);
if any([st1 st2 st3] == 127)
    cleanup_cps96_TI();
    error(['calc_der96_TI: the CPS TI programs (tprep96, tdisp96, t',lower(wavetype),'egn96) are not on the PATH']);
end
Rk.c = zeros(1,0); Rk.u = zeros(1,0); Rk.eig = [];
try
    [~, rk] = read_der96_TI(derfile,wavetype,0);
catch
    return % no valid output at this period (treated as no roots)
end
nr = size(rk.phv,2);
if nr == 0
    return
end
Rk.c = rk.phv(1,:); Rk.u = rk.grv(1,:);
s.iper = ones(1,nr); s.imode = 0:nr-1;
Rk.eig = read_der96_TI(derfile,wavetype,s);
end

function sel = select_branch(Ri, nfreq, T, nmode, dw_frac, U_tol, dup_tol, vs_hs, maxamp, nmax_amp, wavetype, isol)
% Select root nmode (or, with maxamp, the branch among modes 0..nmax_amp
% with the largest surface amplitude) at one period by following branches
% through the computed frequencies (see ROOT VALIDATION).
% Ri = {w, w(1+dw), w(1-dw)}. Returns 1-based root indices jz, jp, jm (NaN
% if none) into Ri{f}, and the selected mode number.
sel = struct('jz',NaN,'jp',NaN,'jm',NaN,'repaired',false,'nspur',0,'mode',NaN);
w = 2*pi/T;
[cz, uz, jz, ndz] = getroots(Ri{1}, dup_tol, vs_hs);
[cp, up, jp]      = getroots(Ri{2}, dup_tol, vs_hs);
if nfreq >= 3
    [cm, um, jm] = getroots(Ri{3}, dup_tol, vs_hs);
else
    cm = zeros(1,0); um = cm; jm = cm;
end

% links between adjacent frequencies
Lzp = linkroots(w, 1, cz, uz, 1+dw_frac, cp, up, U_tol);
Lzm = linkroots(w, 1, cz, uz, 1-dw_frac, cm, um, U_tol);
Lmp = linkroots(w, 1-dw_frac, cm, um, 1+dw_frac, cp, up, U_tol);

% branches through the requested frequency
br = zeros(0,4); % [c, kz, km, kp] (indices into the deduplicated lists, NaN if absent)
usedp = false(size(cp)); usedm = false(size(cm));
for k = 1:length(cz)
    kp = Lzp(k); km = Lzm(k);
    accept = ~isnan(kp) || ~isnan(km) || (isempty(cp) && isempty(cm));
    if accept
        br(end+1,:) = [cz(k), k, km, kp]; %#ok<AGROW>
        if ~isnan(kp); usedp(kp) = true; end
        if ~isnan(km); usedm(km) = true; end
    else
        br(end+1,:) = [cz(k), -k, NaN, NaN]; %#ok<AGROW> % rejected (kept to count)
    end
end
% branches missed at the requested frequency but found at both neighbors
for km = 1:length(cm)
    kp = Lmp(km);
    if ~usedm(km) && ~isnan(kp) && ~usedp(kp)
        br(end+1,:) = [(cm(km)+cp(kp))/2, NaN, km, kp]; %#ok<AGROW>
        usedp(kp) = true;
    end
end
[~, Is] = sort(br(:,1));
br = br(Is,:);
good = find(isnan(br(:,2)) | br(:,2) > 0);
if maxamp
    ncand = min(numel(good), nmax_amp+1);
    if ncand == 0
        sel.nspur = ndz + sum(br(:,2) < 0);
        return
    end
    exc = nan(1,ncand);
    for jc = 1:ncand
        bc = br(good(jc),:);
        if ~isnan(bc(2))
            exc(jc) = surfamp(Ri{1}, jz(bc(2)), wavetype, isol);
        else % missed at w: mean of the neighbors
            exc(jc) = mean([surfamp(Ri{3}, jm(bc(3)), wavetype, isol), ...
                            surfamp(Ri{2}, jp(bc(4)), wavetype, isol)]);
        end
    end
    [~, isel] = max(exc);
    nsel = isel - 1;
else
    if numel(good) < nmode+1
        sel.nspur = ndz + sum(br(:,2) < 0);
        return
    end
    nsel = nmode;
end
sel.mode = nsel;
b = br(good(nsel+1),:);
sel.nspur = sum(br(1:good(nsel+1),2) < 0) + ndz;
if ~isnan(b(2))
    sel.jz = jz(b(2));
else
    sel.repaired = true;
end
if ~isnan(b(3)); sel.jm = jm(b(3)); end
if ~isnan(b(4)); sel.jp = jp(b(4)); end
end

function a = surfamp(Rk, j, wavetype, isol)
% surface amplitude (source and receiver at the top of the solid) of root j
% (1-based CPS root index) in the per-period output Rk:
% AL*ut(top)^2 (Love) or AR*uz(top)^2 (Rayleigh)
if isempty(Rk) || isempty(Rk.eig) || j > length(Rk.eig.amp)
    a = NaN; return
end
e = Rk.eig;
if strcmp(wavetype,'L')
    a = e.amp(j) * e.ut(isol,j)^2;
else
    a = e.amp(j) * e.uz(isol,j)^2;
end
end

function eig = applyQ(eig, model, qa, qb, fref)
eig.phv_elastic = eig.phv;
eig.gamma = zeros(size(eig.phv));
if ~any(qa) && ~any(qb)
    return
end
for ii = 1:length(eig.phv)
    c0 = eig.phv(ii);
    if isnan(c0); eig.gamma(ii) = NaN; continue; end
    f = 1/eig.periods(ii); w = 2*pi*f;
    sumS = eig.dcdvsv(:,ii).*model(:,3) + eig.dcdvsh(:,ii).*model(:,6);
    sumP = eig.dcdvpv(:,ii).*model(:,2) + eig.dcdvph(:,ii).*model(:,5);
    sumS(isnan(sumS)) = 0; sumP(isnan(sumP)) = 0;
    x = sum(sumS.*qb + sumP.*qa);
    eig.gamma(ii) = 0.5*w*x/(c0*c0);
    eig.phv(ii) = c0 + x*log(f/fref)/pi;
end
end

function eig = pickroots(R, f, j, tmpl, n)
% one column per requested period: root j(ii) (1-based) of R{ii,f}
eig = struct();
if isempty(tmpl)
    eig.phv = nan(1,n); eig.grv = nan(1,n);
    return
end
flds = fieldnames(tmpl);
for jj = 1:length(flds)
    if any(strcmp(flds{jj},{'periods','mode'})); continue; end
    eig.(flds{jj}) = nan(size(tmpl.(flds{jj}),1), n);
end
for ii = 1:n
    if isnan(j(ii)) || isempty(R{ii,f}) || isempty(R{ii,f}.eig) || j(ii) > length(R{ii,f}.eig.phv)
        continue
    end
    for jj = 1:length(flds)
        if any(strcmp(flds{jj},{'periods','mode'})); continue; end
        eig.(flds{jj})(:,ii) = R{ii,f}.eig.(flds{jj})(:,j(ii));
    end
end
end

function [c, u, j, ndup] = getroots(Rk, dup_tol, vs_hs)
% roots at one period/frequency, sorted by phase velocity, with duplicates
% and roots at the (flattened) half-space S velocity removed (tdisp96 can
% report non-dispersive roots there, at the upper end of its search range).
% j are the original CPS root indices (1-based). ndup counts the removed
% roots.
c = zeros(1,0); u = c; j = c; ndup = 0;
if isempty(Rk) || isempty(Rk.c)
    return
end
j = find(~isnan(Rk.c));
c = Rk.c(j); u = Rk.u(j);
if vs_hs > 0
    iedge = abs(c/vs_hs - 1) < 1e-4;
    c(iedge) = []; u(iedge) = []; j(iedge) = [];
end
[c, Is] = sort(c); u = u(Is); j = j(Is);
keep = true(size(c));
for k = 2:length(c)
    if abs(c(k)-c(k-1)) <= dup_tol*c(k)
        keep(k) = false;
    end
end
ndup = sum(~keep);
c = c(keep); u = u(keep); j = j(keep);
end

function L = linkroots(w, fa, ca, ua, fb, cb, ub, U_tol)
% For each root a (at angular frequency w*fa) the index of the root b (at
% w*fb) on the same branch: the secant group velocity between them must
% match the mean of their group velocities within U_tol. NaN if there is no
% such root.
L = nan(size(ca));
if isempty(cb)
    return
end
for k = 1:length(ca)
    dk = w*fb./cb - w*fa/ca(k);
    Ufd = w*(fb-fa)./dk;
    e = abs(Ufd - (ua(k)+ub)/2) / abs(ua(k));
    [emin, jb] = min(e);
    if emin <= U_tol
        L(k) = jb;
    end
end
end
