%DEMO_BOUNDARY  Stability boundary of the minimal model at d = 0.0027.
%   Compares the closed form g_crit = 2 w_g / (3 sqrt(Theta_theta)), the adiabatic
%   average of the fixed-amplitude Floquet growth rate, and a short direct
%   estimate of the largest Lyapunov exponent along stochastic theta paths.

c = tg_const();
d = 0.0027;
Theta = c.sigma^2/(2*d*c.wT);
gcrit = 2*c.wG/(3*sqrt(Theta));

rgrid = [linspace(0, 0.89, 90), linspace(0.9, 1.2, 301), linspace(1.21, 20, 940)];
lam2 = tg_lambda2(rgrid, d);                  % growth rate vs theta amplitude / a_c
gAd = tg_adiabatic_gstar(d, rgrid, lam2);
fprintf('closed form g_crit = %.2f, adiabatic Floquet average g* = %.2f\n', gcrit, gAd);

gs = [10 15 20 25 30];
lam = tg_lyapunov_grid(gs, d, 1:10, 10, 50);  % 10 theta paths, 50 s each
fprintf('g = %5.1f   mean largest Lyapunov exponent = %+.2f 1/s\n', [gs; mean(lam, 3).']);
fprintf('(baseline decay -G_gamma0/2 = %.2f 1/s)\n', -d*c.wG/2);
