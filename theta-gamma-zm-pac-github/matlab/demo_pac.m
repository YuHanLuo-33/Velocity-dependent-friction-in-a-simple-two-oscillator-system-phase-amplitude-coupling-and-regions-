%DEMO_PAC  Phase-amplitude coupling of the minimal model against the closed form.
%   Simulates a few records at g_vv = 8.2, d = 0.0027, measures the detrended
%   Tort modulation index of each record and compares the median with
%   MI(kappa*sqrt(ln 2)). Then runs the within-cycle surrogate test on one record.
%   Usage: run from this folder, or addpath it first.

c = tg_const();
g = 8.2; d = 0.0027; nRec = 20;
out = tg_sim(g*ones(1, nRec), d*ones(1, nRec), 1:nRec, 30, 10, ...
    struct('H', 1e12, 'storeMask', true(1, nRec)));

mi = nan(1, nRec); phase = nan(1, nRec);
for i = 1:nRec
    if isfinite(out.cross(i)), continue, end
    P = tg_pac_record(out.pos(:, i, 1), out.pos(:, i, 2), c.h);
    if P.valid
        mi(i) = P.detMI; phase(i) = P.prefPhaseDet;
    end
end

Theta = c.sigma^2/(2*d*c.wT);                 % theta effective temperature
kappa = g*sqrt(2*Theta)/(2*c.wT);             % PAC parameter
fprintf('kappa = %.2f\n', kappa);
fprintf('median MI over %d records: %.3f   closed form MI(kappa sqrt(ln 2)): %.3f\n', ...
    nnz(isfinite(mi)), median(mi, 'omitnan'), tg_mi_vm(kappa*sqrt(log(2))));
fprintf('mean preferred phase: %.3f rad (0 = theta peak)\n', angle(mean(exp(1i*phase), 'omitnan')));

k = find(isfinite(mi), 1);
P = tg_pac_record(out.pos(:, k, 1), out.pos(:, k, 2), c.h);
S = tg_surrogate_test(P, 199, RandStream('mt19937ar', 'Seed', 1));
fprintf('record %d: MI %.3f, surrogate p = %.3f (199 within-cycle shifts)\n', k, P.detMI, S.det.p);
