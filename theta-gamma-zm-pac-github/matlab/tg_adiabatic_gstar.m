function gs = tg_adiabatic_gstar(d, rgrid, lam2)
%TG_ADIABATIC_GSTAR  Coupling at which the Rayleigh average of Lambda_2/2 over theta
%   amplitudes crosses zero (adiabatic stability boundary).
c = tg_const();
sx = sqrt(c.sigma^2/(2*d*c.wT))/c.wT;
a = linspace(0, 12*sx, 4001);
pa = a/sx^2.*exp(-a.^2/(2*sx^2));
gList = 0.5:0.01:60;
vals = zeros(size(gList));
rgrid = rgrid(:); lam2 = lam2(:);
for k = 1:numel(gList)
    rr = a*gList(k)*c.wT/(2*c.wG);
    mu = 0.5*interp1(rgrid, lam2, rr, 'linear', lam2(end));
    vals(k) = trapz(a, pa.*mu);
end
k = find(vals > 0, 1);
if isempty(k) || k == 1
    gs = NaN;
else
    gs = gList(k-1) + (gList(k) - gList(k-1))*(-vals(k-1))/(vals(k) - vals(k-1));
end
end
