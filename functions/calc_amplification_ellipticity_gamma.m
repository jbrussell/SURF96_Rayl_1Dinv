% Calculate Rayleigh wave eigenfunctions 
% 	model: [thickness vp vs]
% 	vec_T: period 
%   nmode: mode branch of interest (0 = fund.)
%   disp_type: 'C'=phase veloc. 'U'=group veloc (OPTIONAL, defaults to C)
%
% jbrussell 8/2023
% 
function disper = calc_amplification_ellipticity_gamma(vec_T,model,wavetype,varargin)

% jbrussell 6/16/2020: Add timeout string to surf96 inversion call to avoid
% hanging. Get gtimeout on mac by installing coreutils with homebrew
% >$ brew install coreutils
% timeoutstr = 'ulimit -t 20; '; % timeout after 20 seconds of CPU time
timeoutstr = 'gtimeout 20 '; % timeout after 20 seconds of "wall" time

system('surf96 39'); % clean up
system('rm start.mod');

if nargin==3
    nmode = 0; % default, fundamental mode
    fref = 1; % default, reference frequency
elseif nargin==4
	nmode = varargin{1};
    fref = 1; % default, reference frequency
elseif nargin==5
    nmode = varargin{1};
    fref = varargin{2}; % reference frequency
end

nmodes = nmode + 1; % number of total modes to calculate

% Make model file
writemod_surf96(model,'start.mod',fref);

% run sprep96 to generate dispersion control file (outputs sdisp96.dat)
fid = fopen('periods.txt','w');
for iper = 1:length(vec_T)
    fprintf(fid,'%f\n',vec_T(iper));
end
fclose(fid);
[~, log] = system([timeoutstr,'sprep96 -M start.mod -NMOD ',num2str(nmodes),' -',wavetype,' -PARR periods.txt']);

% run sdisp96 to generate dispersion binary file (outputs sdisp96.ray or sdisp96.lov)
[~, log] = system([timeoutstr,'sdisp96']);

% calculate eigenfunctions (outputs sregn96.der)
[~, log] = system([timeoutstr,'s',lower(wavetype),'egn96']);

% write dispersion to output file (SREGN.ASC or SLEGN.ASC)
[~, log] = system([timeoutstr,'sdpegn96 -',wavetype,' -C -U -G -A0 -E -ASC']);

disper = read_amp_ellip_gam96(['S',wavetype,'EGN.ASC'],wavetype,nmode);

delete periods.txt start.mod sdisp96.dat sdisp96.ray sdisp96.lov sregn96.der slegn96.der sregn96.egn slegn96.egn SREGNE.PLT SLEGNE.PLT SREGN.ASC SLEGN.ASC SRDER.PLT SLDER.PLT SRDER.TXT SLDER.TXT

end


