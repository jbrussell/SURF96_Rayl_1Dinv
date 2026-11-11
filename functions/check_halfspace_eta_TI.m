% Warn if a TI model has eta ~= 1 in the half-space (Rayleigh waves only)
%
% check_halfspace_eta_TI(model, wavetype)
%
% With eta ~= 1 in the half-space, the CPS program tdisp96 can return a
% spurious, nearly non-dispersive Rayleigh root at low phase velocity
% (e.g. ~1.4 km/s with an unphysical group velocity), which is then reported
% as the fundamental mode. It is triggered when the root search starts below
% that velocity, i.e. for models with a water layer or slow sediments.
% eta ~= 1 in the overlying layers is fine. Setting eta = 1 in the
% half-space avoids the problem.
%
% jbrussell 11/2026
%
function check_halfspace_eta_TI(model, wavetype)

if strcmpi(wavetype,'R') && size(model,2) >= 7 && abs(model(end,7)-1) > 1e-3
    warning('TI:halfspaceEta', ['Half-space eta = %g ~= 1. CPS tdisp96 can return a spurious ', ...
        'low-velocity Rayleigh root for this case (especially with a water layer or slow ', ...
        'sediments). Consider using eta = 1 in the half-space.'], model(end,7));
end

end
