function lam2 = tg_lambda2(r, d)
%TG_LAMBDA2  Growth rate Lambda_2 (1/s) of gamma energy for a theta cycle of fixed
%   amplitude a = r * a_c, a_c = 2 w_g / (g w_th); depends only on r and d.
%   Fundamental matrix over one theta period (RK4), as tg_floquet2.
c = tg_const();
p = struct('fTheta', c.fTheta, 'fGamma', c.fGamma, 'dGamma', d);
g = 1;                                   % any g: R = r a_c gives g w_th R = 2 w_g r
lam2 = tg_floquet2(g, r(:)*2*c.wG/(g*c.wT), p);
end
