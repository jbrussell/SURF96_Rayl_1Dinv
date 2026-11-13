% Common setup for the benchmark_TI_*.m scripts (called by each of them)
%
% Sets the PATH to the CPS binaries, adds ./functions to the MATLAB/Octave
% path, and creates the output directory for figures and the summary log.
%
% The CPS binary directory defaults to ./bin_v3.30/ and can be overridden
% with the environment variable CPS_BIN, e.g. for a Linux build:
%   >$ CPS_BIN=/path/to/PROGRAMS.330/bin octave benchmark_TI_run_all.m
%
% The *_TI functions call the CPS TI programs tprep96, tdisp96, tregn96,
% tlegn96 and tdpegn96 (in addition, the isotropic comparisons use surf96,
% srfker96, sprep96, sdisp96, sregn96, slegn96, sdpegn96 and sdpder96).
% Like the isotropic wrappers they call "gtimeout" (brew install coreutils
% on mac; on Linux, link gtimeout to timeout).
%
% jbrussell 11/2026
%
path2BIN = getenv('CPS_BIN');
if isempty(path2BIN)
    path2BIN = './bin_v3.30/'; % path to CPS binaries
    system('chmod +x ./bin_v3.30/*'); % Make binary files executable
end
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
    setenv('PATH', [path2BIN,':',PATH]);
end
addpath('./functions/')

figdir = './benchmark_TI_figs/';
if ~exist(figdir,'dir')
    mkdir(figdir);
end
logfile = [figdir,'benchmark_TI_summary.txt'];

% Run headless in Octave if there is no display
if exist('OCTAVE_VERSION','builtin') && isempty(getenv('DISPLAY'))
    set(0,'defaultfigurevisible','off');
end
