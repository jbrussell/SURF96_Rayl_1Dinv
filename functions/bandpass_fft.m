% Zero-phase band-pass filter in the frequency domain (cosine tapers)
%
%   y = bandpass_fft(x, dt, Tmin, Tmax)
%   y = bandpass_fft(x, dt, Tmin, Tmax, taper)
%
% x:     time series [npts x ntrace] (filtered along columns)
% dt:    sampling interval (s)
% Tmin, Tmax: pass band in period (s)
% taper: width of the cosine tapers as a fraction of the corner
%        frequencies (default 0.2)
%
% The trace mean is removed and the series zero-padded to twice its length
% before filtering.
%
% jbrussell 11/2026
%
function y = bandpass_fft(x, dt, Tmin, Tmax, taper)
if nargin < 5; taper = 0.2; end
if isvector(x); x = x(:); end
[n, ntr] = size(x);
N = 2^nextpow2(2*n);
f1 = 1/Tmax; f2 = 1/Tmin;
f = (0:N-1)'/(N*dt);
f(f > 1/(2*dt)) = f(f > 1/(2*dt)) - 1/dt;
af = abs(f);
w = zeros(N,1);
w(af >= f1 & af <= f2) = 1;
df1 = taper*f1; df2 = taper*f2;
i = af > f1-df1 & af < f1;
w(i) = 0.5*(1 - cos(pi*(af(i) - (f1-df1))/df1));
i = af > f2 & af < f2+df2;
w(i) = 0.5*(1 + cos(pi*(af(i) - f2)/df2));
X = fft(x - repmat(mean(x,1), n, 1), N);
y = real(ifft(X .* repmat(w, 1, ntr)));
y = y(1:n,:);
end
