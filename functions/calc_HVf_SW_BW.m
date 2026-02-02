function [freq,hvsr,vel] = calc_HVf_SW_BW(model, nf,fmin,fmax,nmr,nml,nks)
% Calculate HVSR considering surface waves and body waves
%
% Usage: HV [OPTIONS]
% 
%  OPTIONS FOR CONTROLING HV:
% 
%  -fmin X : Minimum frequency
%  -fmax X : Maximum frequency
%  -nf N   : Number of frequencies
%  -nmr  N : Max. of Rayleigh modes to be considered
%  -nml  N : Max. of Love modes to be considered
%  -prec X : Rel. precission in slowness (default: 1E-4 per cent)
%  -nks  N : Number of k values for numeric integrals
%            Use 0 or missing flag to skip BW calculation
%  -apsv X : Attenuation to stabilize PSV. w -> w-I*apsv*w (DEFAULT = 0)
%  -ash  X : Attenuation to stabilize SH.  w -> w-I*ash*w  (DEFAULT = 0)
%  -logsam  : Regular sampling in log[frequency] (default is regular in frequency)
%  -ff <file> frequency list (sorted) from file
% 
%  INPUT MODEL PROPERTIES:
% 
%  -nlyr N : Number of layer (including halfspace)
%  -dens Density
%  -thk  Layer thicknesses
%  -vp   P wave velocity
%  -vs   S wave velocity
%  -f <file> Model from file (see below for format)
% 
%  OUTPUTS SELECTION:
% 
%  -hv     : Outputs ( freq., H/V) pairs -> HV.dat
%  -ph     : Outputs phase slowness -> Rph.dat,Lph.dat containing:
%            Number of frequencies (N), number of modes (NM), list of N x NM
%            slowness values for increasing frequencies and mode number, list of
%            N x NM logical values. False (F) and 0 slowness indicate that the
%            mode does not exist at that frequency or could not be calculated.
%            The inner loop runs on the frequency.
%  -gr     : Outputs group velocities -> Rgr.dat,Lgr.dat. Same file structure.
%  -rep    : Writes a report if computation fails and details of calculations ->
%          -> report.dat,GBW.dat (body waves),G123.dat (surface waves),ERR.DAT 
%  The last flag listed of these four will write into the std output instead of
%  using the specified file
% 
%  FORMAT FOR LAYERED MODEL FILES :
% 
%  Line 1    <number of layers including halfspace>
%  Line 2    <thickness(m)> <Vp(m / s)> <Vs(m / s)> <Density(kg / m3)>
%  ...
%  Line N    0 <Vp(m / s)> <Vs(m / s)> <Density(kg / m3)>
% 
%
% USAGE(examples)
% 
% HV.exe -f model.txt -nf 100 -fmin 0.1 -fmax 10 -logsam -nmr 20 -nml 20 -ph -hv  > HV.dat
%

fid = fopen('model.txt','w');
nlyr = size(model,1);
thk = model(:,1)*1000;
vp = model(:,2)*1000;
vs = model(:,3)*1000;
dens = model(:,4)*1000;
fprintf(fid,'%d\n',nlyr);
for ii = 1:nlyr
    fprintf(fid,'%8.2f %8.2f %8.2f %8.2f\n',thk(ii),vp(ii),vs(ii),dens(ii));
end
fclose(fid);

[~,log] = system(['HVf -f model.txt -nf ',num2str(nf),' -fmin ',num2str(fmin),' -fmax ',num2str(fmax),' -nmr ',num2str(nmr),' -nml ',num2str(nml),' -nks ',num2str(nks),' -logsam -ph -gr -hv  > HV.dat']);

% Load HVSR
dat = load('HV.dat');
dat = reshape(dat,2,length(dat)/2)';
freq = dat(:,1);
hvsr = dat(:,2);

% Initialize Rayleigh and Love Group Velocities
Rph=[];
Lph=[];
Rgr=[];
Lgr=[];

% Load Rayleigh and Love Phase velocity
if exist('Rph.dat')
    [Rph] = load_v('Rph.dat');
end
if exist('Lph.dat')
    [Lph] = load_v('Lph.dat');
end

% Load Rayleigh and Love Group velocity
if exist('Rgr.dat')
    [Rgr] = load_v('Rgr.dat');
end
if exist('Lgr.dat')
    [Lgr] = load_v('Lgr.dat');
end

vel.Rph = Rph;
vel.Lph = Lph;
vel.Rgr = Rgr;
vel.Lgr = Lgr;


delete HV.dat Lgr.dat Lph.dat model.txt Rgr.dat Rph.dat

end

% Load Velocity
function [v] = load_v(fname)
    fid = fopen(fname);
    row = fgetl(fid);
    vals = str2num(row);
    Nf = vals(1);
    Nmode = vals(2);
    
    row = fgetl(fid);
    s = str2num(row);
    s = reshape(s,Nf,Nmode);
    v = 1./s / 1000; % km/s
    v(isinf(v)) = nan;
    
    fclose(fid);
end




