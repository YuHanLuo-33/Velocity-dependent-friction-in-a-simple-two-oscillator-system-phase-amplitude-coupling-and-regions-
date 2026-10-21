function P = tg_pac_record(xTheta, xGamma, h)
%TG_PAC_RECORD  PAC of one record (see the paper's appendix).
%   Theta phase and gamma envelope from FFT analytic signals of the centred
%   positions; the envelope is divided by one fitted log-linear trend (fit over
%   samples above 1e-3 of the maximum); Tort MI with 18 phase bins.
%   P.valid is false for records with non-finite values or empty bins.
if nargin < 3, c = tg_const(); h = c.h; end
K = 18;
P = struct('valid', false, 'rawMI', NaN, 'detMI', NaN, 'prefPhaseRaw', NaN, ...
    'prefPhaseDet', NaN, 'trendFactor', NaN, 'trendSlope', NaN);
th = double(xTheta(:));  gm = double(xGamma(:));
if any(~isfinite(th)) || any(~isfinite(gm)), return, end
N = numel(th);
t = h*(1:N).';
phase = angle(tg_analytic(th - mean(th)));
amp = abs(tg_analytic(gm - mean(gm)));
mask = isfinite(amp) & amp > 1e-3*max(amp);
if nnz(mask) < 10, return, end
cf = polyfit(t(mask), log(amp(mask)), 1);
trend = max(exp(cf(2) + cf(1)*t), eps);
det = amp./trend;
bins = tg_bin_index(phase, K);
counts = accumarray(bins, 1, [K 1]);
if any(counts == 0), return, end
mraw = accumarray(bins, amp, [K 1])./counts;
mdet = accumarray(bins, det, [K 1])./counts;
if any(mraw <= 0) || any(mdet <= 0), return, end
edges = linspace(-pi, pi, K + 1);
centers = edges(1:end-1).' + pi/K;
P.valid = true;
P.phase = phase; P.bins = bins; P.counts = counts; P.amp = amp; P.det = det;
P.rawMI = tg_mi_from_means(mraw);
P.detMI = tg_mi_from_means(mdet);
P.prefPhaseRaw = angle(sum(mraw.*exp(1i*centers)));
P.prefPhaseDet = angle(sum(mdet.*exp(1i*centers)));
P.trendFactor = exp(abs(cf(1))*(t(end) - t(1)));
P.trendSlope = cf(1);
end
