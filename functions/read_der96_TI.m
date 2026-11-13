% Read the binary eigenfunction/partial-derivative file written by the CPS
% transversely isotropic programs "tregn96 -DER" (tregn96.der) or
% "tlegn96 -DER" (tlegn96.der).
%
% This is the TI analog of read_eig96.m. Rather than parsing the text file
% written by tdpder96 -TXT (which only carries 3 significant figures), it
% reads the Fortran unformatted (binary) file directly, so eigenfunctions
% and kernels are returned at full single precision. This extra precision
% is needed to compute group-velocity kernels by finite differences
% (see calc_kernel96_TI.m).
%
%   eig = read_der96_TI(filename, wavetype)
%   eig = read_der96_TI(filename, wavetype, nmode)
%   [eig, raw] = read_der96_TI(filename, wavetype, sel)
%
%   filename: 'tregn96.der' or 'tlegn96.der'
%   wavetype: 'R' or 'L'
%   nmode:    mode branch of interest (0 = fund.) (OPTIONAL, default 0)
%   sel:      (alternative to nmode) structure selecting one root per output
%             column: sel.iper = indices into raw.periods, sel.imode = CPS
%             root index (0 = first root found; NaN = none). Used by
%             calc_der96_TI.m, which discards spurious roots.
%
% Output structure (each matrix is [nlayer x nperiod]); periods at which the
% requested mode does not exist are filled with NaN:
%   eig.periods, eig.mode, eig.phv, eig.grv, eig.gamma, eig.amp (AR or AL),
%   eig.ellip (Rayleigh ellipticity UR/UZ at the surface, as tdpegn96 -E)
%   Rayleigh: eig.ur, eig.tr, eig.uz, eig.tz
%   Love:     eig.ut, eig.tt
%   Phase-velocity partial derivatives (layer kernels, eta held fixed):
%     eig.dcdh, eig.dcdvsv, eig.dcdvsh, eig.dcdrho
%     Rayleigh only: eig.dcdvpv, eig.dcdvph, eig.dcdeta
%     (For Love waves dcdvpv, dcdvph and dcdeta are returned as zeros.)
%
% raw: all roots found by CPS: raw.periods (1 x nper), raw.phv and raw.grv
%      ([nper x maxroots], NaN padded)
%
% The binary layout follows tio.f (getmdt / gethed / getdrt) in CPS
% VOLIII/src. The isotropic files written by "sregn96 -DER" / "slegn96 -DER"
% (sio.f) have the same layout and can also be read: dc/dVp and dc/dVs are
% then returned as dcdvpv and dcdvsv (dcdvph, dcdvsh and dcdeta are NaN). Fortran sequential unformatted records with 4-byte
% little-endian record markers (gfortran default) are assumed.
%
% jbrussell 11/2026
%
function [eig, raw] = read_der96_TI(filename,wavetype,varargin)

if nargin==2
    nmode = 0; % default, fundamental mode
elseif nargin==3
    nmode = varargin{1};
end
wavetype = upper(wavetype);

fid = fopen(filename,'r','ieee-le');
if fid < 0
    error([filename,' does not exist']);
end
recs = {};
while true
    n1 = fread(fid,1,'int32');
    if isempty(n1)
        break
    end
    bytes = fread(fid,n1,'uint8=>uint8');
    n2 = fread(fid,1,'int32');
    if isempty(n2) || n2~=n1 || numel(bytes)~=n1
        break % truncated or corrupt record
    end
    recs{end+1} = bytes(:)'; %#ok<AGROW>
end
fclose(fid);

if numel(recs) < 5
    error([filename,' does not appear to be a valid TI eigenfunction file']);
end

% ---- File header (getmdt) ----
ipar = double(typecast(recs{2},'int32'));
mmax = double(typecast(recs{4}(1:4),'int32'));
irec = 6;

% ---- Loop over periods ----
periods = [];
blocks = {};
while irec <= numel(recs)
    r = recs{irec};
    irec = irec + 1;
    if numel(r) < 12
        break
    end
    ifunc = double(typecast(r(1:4),'int32'));
    kmode = double(typecast(r(5:8),'int32'));
    t0 = double(typecast(r(9:12),'single'));
    if ifunc <= 0
        break
    end
    blk = {};
    for j = 1:kmode
        b = struct();
        hvals = double(typecast(recs{irec},'single'));
        irec = irec + 1;
        b.wvno = hvals(1);
        b.u = hvals(2);
        b.gamma = hvals(3);
        b.amp = hvals(end); % AR (Rayleigh) or AL (Love) energy factor
        if strcmp(wavetype,'R')
            b.ellip = hvals(14); % ellipticity (UR/UZ) at the free surface
        else
            b.ellip = NaN;
        end
        if strcmp(wavetype,'L')
            if ipf(ipar,4)
                [b.ur, irec] = getarr(recs, irec); [b.tr, irec] = getarr(recs, irec);
            end
            if ipf(ipar,5);  [b.dcdh, irec]  = getarr(recs, irec); end
            if ipf(ipar,7);  [b.dcdbv, irec] = getarr(recs, irec); end
            if ipf(ipar,8);  [b.dcdr, irec]  = getarr(recs, irec); end
            if ipf(ipar,11); [b.dcdbh, irec] = getarr(recs, irec); end
        else
            if ipf(ipar,4)
                [b.ur, irec] = getarr(recs, irec); [b.tr, irec] = getarr(recs, irec);
                [b.uz, irec] = getarr(recs, irec); [b.tz, irec] = getarr(recs, irec);
            end
            if ipf(ipar,5);  [b.dcdh, irec]  = getarr(recs, irec); end
            if ipf(ipar,6);  [b.dcdav, irec] = getarr(recs, irec); end
            if ipf(ipar,7);  [b.dcdbv, irec] = getarr(recs, irec); end
            if ipf(ipar,8);  [b.dcdr, irec]  = getarr(recs, irec); end
            if ipf(ipar,9);  [b.dcdah, irec] = getarr(recs, irec); end
            if ipf(ipar,10); [b.dcdn, irec]  = getarr(recs, irec); end
            if ipf(ipar,11); [b.dcdbh, irec] = getarr(recs, irec); end
        end
        blk{j} = b; %#ok<AGROW>
    end
    periods(end+1) = t0; %#ok<AGROW>
    blocks{end+1} = blk; %#ok<AGROW>
end
nper = numel(periods);

% ---- Table of all roots ----
maxroots = max([0, cellfun(@numel,blocks)]);
raw.periods = periods(:)';
raw.phv = nan(nper,maxroots);
raw.grv = nan(nper,maxroots);
for ip = 1:nper
    for j = 1:numel(blocks{ip})
        raw.phv(ip,j) = 2*pi/periods(ip) / blocks{ip}{j}.wvno;
        raw.grv(ip,j) = blocks{ip}{j}.u;
    end
end

% ---- Which root goes in each output column ----
if isstruct(nmode)
    iper = nmode.iper(:)';
    imode = nmode.imode(:)';
    eig.mode = NaN;
else
    iper = 1:nper;
    imode = nmode*ones(1,nper);
    eig.mode = nmode;
end
nout = numel(iper);

% ---- Arrange output ----
eig.periods = periods(iper);
eig.phv = nan(1,nout);
eig.grv = nan(1,nout);
eig.gamma = nan(1,nout);
eig.amp = nan(1,nout);
eig.ellip = nan(1,nout);
if strcmp(wavetype,'L')
    flds_in  = {'ur','tr','dcdh','dcdbv','dcdbh','dcdr'};
    flds_out = {'ut','tt','dcdh','dcdvsv','dcdvsh','dcdrho'};
else
    flds_in  = {'ur','tr','uz','tz','dcdh','dcdav','dcdah','dcdn','dcdbv','dcdbh','dcdr'};
    flds_out = {'ur','tr','uz','tz','dcdh','dcdvpv','dcdvph','dcdeta','dcdvsv','dcdvsh','dcdrho'};
end
for ii = 1:numel(flds_out)
    eig.(flds_out{ii}) = nan(mmax,nout);
end
for io = 1:nout
    ip = iper(io);
    j = imode(io) + 1;
    if isnan(j) || j > numel(blocks{ip})
        continue % requested mode does not exist at this period
    end
    b = blocks{ip}{j};
    eig.phv(io) = 2*pi/periods(ip) / b.wvno;
    eig.grv(io) = b.u;
    eig.gamma(io) = b.gamma;
    eig.amp(io) = b.amp;
    eig.ellip(io) = b.ellip;
    for ii = 1:numel(flds_in)
        if isfield(b,flds_in{ii})
            eig.(flds_out{ii})(:,io) = b.(flds_in{ii})(:);
        end
    end
end
if strcmp(wavetype,'L')
    % Love waves have no sensitivity to these parameters
    eig.dcdvpv = zeros(mmax,nout);
    eig.dcdvph = zeros(mmax,nout);
    eig.dcdeta = zeros(mmax,nout);
    eig.dcdvpv(:,isnan(eig.phv)) = NaN;
    eig.dcdvph(:,isnan(eig.phv)) = NaN;
    eig.dcdeta(:,isnan(eig.phv)) = NaN;
end

end

function tf = ipf(ipar, k)
% true if output flag k is set (the isotropic s?egn96 files have fewer flags)
tf = numel(ipar) >= k && ipar(k) == 1;
end

function [a, irec] = getarr(recs, irec)
% read the next record as a column vector of single-precision values
a = double(typecast(recs{irec},'single'));
a = a(:);
irec = irec + 1;
end
