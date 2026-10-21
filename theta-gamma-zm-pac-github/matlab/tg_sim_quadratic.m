function [posT, posG, ok] = tg_sim_quadratic(c2, c1, bxx, bvx, d, seeds, tBurn, tRec, H)
%TG_SIM_QUADRATIC  Gamma driven by the general quadratic theta-gamma terms.
%   dv_g = [-(G0 - c2 v_th - c1 x_th) v_g - (w_g^2 + bxx x_th + bvx v_th) x_g] dt + sigma dW_g
%   c2: v_th v_g (g_vv of the paper); c1: x_th v_g; bxx: x_th x_g; bvx: v_th x_g.
%   Theta is unchanged. Same trapezoidal scheme, streams, stationary start and
%   burn-in as tg_sim; theta-dependent coefficients use the step average of theta.
%   Records are not stopped; ok = max(|x_g|, |v_g|) stayed below H (default 1e12).
%   posT, posG: nRec x R (double).
if nargin < 7, tBurn = 30; end
if nargin < 8, tRec = 10; end
if nargin < 9, H = 1e12; end
c = tg_const(); h = c.h;
seeds = seeds(:).';
R = numel(seeds);
[z0, streams, inv] = tg_streams(seeds);
GT = d*c.wT; GG = d*c.wG;
sdVT = c.sigma/sqrt(2*GT); sdVG = c.sigma/sqrt(2*GG);
xT = sdVT/c.wT*z0(1, :); xG = sdVG/c.wG*z0(2, :);
vT = sdVT*z0(3, :);      vG = sdVG*z0(4, :);
qT = 1/(1 + 0.5*h*GT); aT = qT*(1 - 0.5*h*GT); bT = qT*h*c.wT^2; nsT = qT*c.sigma*sqrt(h);
ns = c.sigma*sqrt(h);
nb = round(tBurn/h); nr = round(tRec/h);
posT = zeros(nr, R); posG = zeros(nr, R); mx = zeros(1, R);
step = 0;
while step < nb + nr
    nC = min(1000, nb + nr - step);
    nz = tg_noise_chunk(streams, inv, nC);
    for k = 1:nC
        step = step + 1;
        vTn = aT*vT - bT*xT + nsT*nz(1, :, k);
        xTn = xT + h*vTn;
        cc = GG - 0.5*c2*(vT + vTn) - 0.5*c1*(xT + xTn);
        w2 = c.wG^2 + 0.5*bxx*(xT + xTn) + 0.5*bvx*(vT + vTn);
        vG = ((1 - 0.5*h*cc).*vG - h*w2.*xG + ns*nz(2, :, k))./(1 + 0.5*h*cc);
        xG = xG + h*vG;
        vT = vTn; xT = xTn;
        mx = max(mx, max(abs(xG), abs(vG)));
        j = step - nb;
        if j >= 1
            posT(j, :) = xT; posG(j, :) = xG;
        end
    end
end
ok = mx < H;
end
