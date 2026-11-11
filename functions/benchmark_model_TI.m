% Isotropic oceanic model used by the benchmark_TI_*.m scripts
%
% model = benchmark_model_TI()        -> [H Vp Vs rho]
% model = benchmark_model_TI(1)       -> [H Vp Vs rho Qp Qs]
%
% 3 km water, 0.5 km sediment, 6.5 km crust (1 km layers), mantle to 150 km
% (5 km layers) over a half-space.
%
% jbrussell 11/2026
%
function model = benchmark_model_TI(varargin)

withQ = 0;
if nargin >= 1
    withQ = varargin{1};
end

z = [0 3 3.5 4:1:10 15:5:150]';
dz = [diff(z); 0];
vs = 3.4 + 0.04*(z-4);   vs(z>=15) = 4.4 + 0.002*(z(z>=15)-15);
vp = 6.2 + 0.05*(z-4);   vp(z>=15) = 8.0 + 0.002*(z(z>=15)-15);
rho = 2.8*ones(size(z)); rho(z>=15) = 3.3;
vs(1) = 0;   vp(1) = 1.5; rho(1) = 1.03; % water
vs(2) = 0.8; vp(2) = 2.0; rho(2) = 1.8;  % sediment
model = [dz vp vs rho];

if withQ
    qp = 1000*ones(size(z)); qs = 400*ones(size(z));
    qp(1) = 10000; qs(1) = 0;            % water
    qp(2) = 200; qs(2) = 100;            % sediment
    qp(z>=15 & z<80) = 400; qs(z>=15 & z<80) = 150;   % lid
    qp(z>=80) = 200; qs(z>=80) = 80;     % asthenosphere
    model = [model qp qs];
end

end
