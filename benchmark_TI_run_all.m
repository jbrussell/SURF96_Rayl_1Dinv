% Run all isotropic vs transversely isotropic (TI) benchmarks and print a
% summary. Figures and the summary are written to ./benchmark_TI_figs/
%
%   benchmark_TI_dispersion        dispR/dispL_surf96 vs *_TI
%   benchmark_TI_kernels           calc_kernel96 vs calc_kernel96_TI, and
%                                  finite-difference checks of the TI kernels
%   benchmark_TI_eigenfunctions    calc_eigenfunctions96, amplification/
%                                  ellipticity/gamma, admittance vs *_TI
%   benchmark_TI_anelastic_kernels calc_anelastic_kernel96 vs *_TI, and
%                                  self-consistency/FD checks
%
% To use CPS binaries other than ./bin_v3.30/ set the environment variable
% CPS_BIN before starting MATLAB/Octave (see benchmark_TI_setup.m).
%
% jbrussell 11/2026
%
clear; close all;
if exist('./benchmark_TI_figs/benchmark_TI_summary.txt','file')
    delete('./benchmark_TI_figs/benchmark_TI_summary.txt');
end
% (each script starts with "clear", so they are called one after another)
benchmark_TI_dispersion;
benchmark_TI_kernels;
benchmark_TI_eigenfunctions;
benchmark_TI_anelastic_kernels;
summary = fileread('./benchmark_TI_figs/benchmark_TI_summary.txt');
fprintf('\n\n================ SUMMARY ================\n%s\n', summary);
nfail = length(strfind(summary,'[FAIL]'));
npass = length(strfind(summary,'[PASS]'));
fprintf('%d passed, %d failed\n', npass, nfail);
