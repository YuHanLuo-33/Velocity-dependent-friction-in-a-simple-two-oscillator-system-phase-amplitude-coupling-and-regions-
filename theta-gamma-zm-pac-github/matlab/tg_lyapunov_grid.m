function [lam, flag] = tg_lyapunov_grid(gs, ds, seeds, tSkip, tMeas)
%TG_LYAPUNOV_GRID  Largest Lyapunov exponent of the noise-free gamma equation
%   along stationary stochastic theta paths, for every (g, d, seed).
%   Theta path for a seed is the same as in tg_sim (same streams). The gamma
%   perturbation uses the trapezoidal update and the energy norm
%   sqrt((w_g dx)^2 + dv^2), renormalised every 10 steps; tSkip seconds are
%   discarded and the exponent is averaged over tMeas seconds. Steps with
%   1 + h c/2 <= 0.05 (step too coarse, deep in the unstable region) are flagged
%   and the exponent is set to 50 s^-1.
%   Output lam, flag: nG x nD x nSeed.
if nargin < 4, tSkip = 10; end
if nargin < 5, tMeas = 100; end
c = tg_const();
h = c.h;
LAMCAP = 50; RENORM = 10;
gs = gs(:); ds = ds(:).'; seeds = seeds(:).';
nG = numel(gs); nD = numel(ds); nS = numel(seeds);
D = repelem(ds, nS);                    % record order: seeds fastest within d
S = repmat(seeds, 1, nD);
[z0, streams, inv] = tg_streams(S);
sdVT = c.sigma./sqrt(2*D*c.wT);
xT = sdVT/c.wT.*z0(1, :);
vT = sdVT.*z0(3, :);
GT = D*c.wT;  GG = D*c.wG;
qT = 1./(1 + 0.5*h*GT);
aT = qT.*(1 - 0.5*h*GT);
bT = qT*h*c.wT^2;
nsT = qT*c.sigma*sqrt(h);
hw2 = h*c.wG^2;
dx = ones(nG, numel(D))/(sqrt(2)*c.wG);
dv = ones(nG, numel(D))/sqrt(2);
logs = zeros(nG, numel(D));
flag = false(nG, numel(D));
nSkip = round(tSkip/h);
nTot = nSkip + round(tMeas/h);
g2 = 0.5*gs;
step = 0;
while step < nTot
    nC = min(2000, nTot - step);
    nz = tg_noise_chunk(streams, inv, nC);
    for j = 1:nC
        step = step + 1;
        vTn = aT.*vT - bT.*xT + nsT.*nz(1, :, j);
        cc = GG - g2*(vT + vTn);                 % nG x nRec
        hc = 0.5*h*cc;
        den = 1 + hc;
        flag = flag | den <= 0.05;
        den(den <= 0.05) = 1;
        dv = ((1 - hc).*dv - hw2*dx)./den;
        dx = dx + h*dv;
        vT = vTn;
        xT = xT + h*vT;
        if mod(step, RENORM) == 0
            nr = sqrt((c.wG*dx).^2 + dv.^2);
            nr(~isfinite(nr) | nr <= 0) = 1;
            if step > nSkip
                logs = logs + log(nr);
            end
            dx = dx./nr;
            dv = dv./nr;
        end
    end
end
lam = logs/tMeas;
lam = min(lam, LAMCAP);
lam(flag) = LAMCAP;
lam = reshape(lam, nG, nS, nD);
lam = permute(lam, [1 3 2]);             % nG x nD x nSeed
flag = permute(reshape(flag, nG, nS, nD), [1 3 2]);
end
