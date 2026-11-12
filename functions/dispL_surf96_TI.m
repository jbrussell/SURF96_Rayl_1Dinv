% function of calculating Love wave dispersion for a transversely
% isotropic (TI) model (TI analog of dispL_surf96.m)
%
% 	model: TI model [H Vpv Vsv rho Vph Vsh eta]  (elastic)
%          or [H Vpv Vsv rho Vph Vsh eta Qp Qs]   (Q-corrected, as in dispL_surf96)
%          See writemod_surf96_TI.m; use convert_iso2TI.m to convert an
%          isotropic [H Vp Vs rho (Qp Qs)] model.
% 	vec_T: period
% 	data: [phv];
%   nmode: 0=fund. 1=1st, ... etc. (OPTIONAL, defaults to fund.)
%   disp_type: 'C'=phase veloc. 'U'=group veloc (OPTIONAL, defaults to C)
%   fref: reference frequency for Q (OPTIONAL, defaults to 1 Hz)
%   opts: (OPTIONAL) options passed to calc_der96_TI.m: opts.nmodes =
%         number of modes computed by tdisp96 (default nmode+1), and
%         opts.c_guess = expected phase velocities (faster; see calc_der96_TI.m)
%
% Periods at which the requested mode does not exist are returned as NaN.
%
% Notes:
%  - surf96 does not support TI models, so (unlike dispL_surf96) the
%    elastic case is computed with tprep96/tdisp96/tlegn96 using
%    Q^-1 = 0 in every layer.
%  - If Qp and Qs are provided, the phase velocity includes the causal
%    (physical dispersion) correction computed by tlegn96. As in CPS, the
%    group velocity is not Q-corrected.
%
% jbrussell 11/2026
%
function phv = dispL_surf96_TI(vec_T,model,varargin)

    if nargin==2
        nmode = 0; % default, fundamental mode
        disp_type = 'C'; % default, phase velocity
        fref = 1; % default, reference frequency
    elseif nargin==3
        nmode = varargin{1};
        disp_type = 'C'; % default, phase velocity
        fref = 1; % default, reference frequency
    elseif nargin==4
        nmode = varargin{1};
        disp_type = upper(varargin{2});
        if ~(strcmp(disp_type,'C') || strcmp(disp_type,'U'))
            error('disp_type must be C or U');
        end
        fref = 1; % default, reference frequency
    elseif nargin==5
        nmode = varargin{1};
        disp_type = upper(varargin{2});
        if ~(strcmp(disp_type,'C') || strcmp(disp_type,'U'))
            error('disp_type must be C or U');
        end
        fref = varargin{3};
    elseif nargin==6
        nmode = varargin{1};
        disp_type = upper(varargin{2});
        if ~(strcmp(disp_type,'C') || strcmp(disp_type,'U'))
            error('disp_type must be C or U');
        end
        fref = varargin{3};
    end
    opts = struct();
    if nargin >= 6 && ~isempty(varargin{4})
        opts = varargin{4};
    end

vec_T = vec_T(:);
n_requested = length(vec_T);

if size(model,2)==7
    % No Q provided: elastic calculation (Q^-1 = 0 in all layers)
    model = [model, zeros(size(model,1),2)];
elseif size(model,2)~=9
    error('TI model must have either 7 or 9 columns: [H Vpv Vsv rho Vph Vsh eta (Qp Qs)]');
end

phv = NaN(n_requested, 1);
try
    % Dispersion is read from the binary eigenfunction file, which carries
    % full single precision (the tdpegn96 ASCII table has only 4 decimals)
    disper = calc_eigenfunctions96_TI(vec_T,model,'L',nmode,fref,opts);
catch
    disper = [];
end

if isempty(disper) || ~isfield(disper,'periods') || isempty(disper.periods)
    return
end

T_returned = disper.periods(:);
if strcmp(disp_type,'C')
    v_returned = disper.phv(:);
elseif strcmp(disp_type,'U')
    v_returned = disper.grv(:);
end

% Match returned periods to the requested periods
for ii = 1:n_requested
    [dT, I] = min(abs(T_returned - vec_T(ii)));
    if ~isempty(I) && dT <= max(0.01, 1e-3*vec_T(ii))
        phv(ii) = v_returned(I);
    end
end

phv = phv(:);

end
