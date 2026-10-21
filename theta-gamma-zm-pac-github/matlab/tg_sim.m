function out = tg_sim(g, d, seeds, tBurn, tRec, opt)
%TG_SIM  Simulate many theta-gamma records with the protocol of the paper.
%
%   out = tg_sim(g, d, seeds, tBurn, tRec, opt)
%     g, d, seeds   vectors with one entry per record (records that share a seed
%                   share the initial state and the noise path)
%     tBurn, tRec   burn-in with coupling and analysed record (s)
%     opt.H         threshold on max(|x_theta|,|v_theta|,|x_gamma|,|v_gamma|); a record
%                   that reaches it (or has 1 + h c/2 <= 0.05) is stopped   [1e6]
%     opt.storeMask logical vector: records whose positions are stored      [none]
%     opt.energyEvery  store gamma energy every n record steps (0 = no)     [0]
%
%   Update (semi-implicit; coupling inside the trapezoidal friction step):
%     v_th' = aT v_th - bT x_th + nsT xi_th
%     c     = G_gamma0 - g (v_th + v_th')/2
%     v_g'  = [(1 - h c/2) v_g - h w_g^2 x_g + sigma sqrt(h) xi_g] / (1 + h c/2)
%     x' = x + h v'
%
%   Output: cross (step from simulation start at which the record was stopped,
%   NaN if never), nBurn, nRec, pos (nRec x nStore x 2 single: theta, gamma; NaN
%   after a stop), storeIdx, energy (nEn x R single), emax (1 x R), alive.

if nargin < 6, opt = struct(); end
if ~isfield(opt, 'H'), opt.H = 1e6; end
if ~isfield(opt, 'storeMask'), opt.storeMask = []; end
if ~isfield(opt, 'energyEvery'), opt.energyEvery = 0; end
c = tg_const();
h = c.h;
g = g(:).'; d = d(:).'; seeds = seeds(:).';
R = numel(g);
nBurn = round(tBurn/h);
nRec = round(tRec/h);
nTot = nBurn + nRec;

[z0, streams, inv] = tg_streams(seeds);
GT = d*c.wT;  GG = d*c.wG;
sdVT = c.sigma./sqrt(2*GT);  sdVG = c.sigma./sqrt(2*GG);
xT = sdVT/c.wT.*z0(1, :);  xG = sdVG/c.wG.*z0(2, :);
vT = sdVT.*z0(3, :);       vG = sdVG.*z0(4, :);
qT = 1./(1 + 0.5*h*GT);
aT = qT.*(1 - 0.5*h*GT);
bT = qT*h*c.wT^2;
nsT = qT*c.sigma*sqrt(h);
ns = c.sigma*sqrt(h);
hw2 = h*c.wG^2;

alive = true(1, R);
cross = nan(1, R);
storeMask = false(1, R);
if ~isempty(opt.storeMask), storeMask = logical(opt.storeMask(:).'); end
sidx = find(storeMask);
pos = nan(nRec, numel(sidx), 2, 'single');
ee = opt.energyEvery;
if ee > 0
    energy = nan(floor(nRec/ee), R, 'single');
else
    energy = [];
end
emax = zeros(1, R);

chunk = 1000;
step = 0;
while step < nTot
    nC = min(chunk, nTot - step);
    nz = tg_noise_chunk(streams, inv, nC);
    for j = 1:nC
        step = step + 1;
        vTn = aT.*vT - bT.*xT + nsT.*nz(1, :, j);
        cc = GG - 0.5*g.*(vT + vTn);
        den = 1 + 0.5*h*cc;
        vG = ((1 - 0.5*h*cc).*vG - hw2*xG + ns*nz(2, :, j))./den;
        xG = xG + h*vG;
        vT = vTn;
        xT = xT + h*vT;
        mx = max(max(abs(vG), abs(xG)), max(abs(vT), abs(xT)));
        bad = alive & (mx >= opt.H | den <= 0.05 | ~isfinite(mx));
        if any(bad)
            cross(bad) = step;
            alive(bad) = false;
            xG(bad) = 0;
            vG(bad) = 0;
        end
        if step > nBurn
            r = step - nBurn;
            if ~isempty(sidx)
                pos(r, :, 1) = xT(sidx);
                gs = xG(sidx);
                gs(~alive(sidx)) = NaN;
                pos(r, :, 2) = gs;
            end
            E = 0.5*(vG.^2 + (c.wG*xG).^2);
            emax(alive) = max(emax(alive), E(alive));
            if ee > 0 && mod(r, ee) == 0
                E(~alive) = NaN;
                energy(r/ee, :) = E;
            end
        end
    end
end
out = struct('cross', cross, 'nBurn', nBurn, 'nRec', nRec, 'pos', pos, ...
    'storeIdx', sidx, 'energy', energy, 'emax', emax, 'alive', alive);
end
