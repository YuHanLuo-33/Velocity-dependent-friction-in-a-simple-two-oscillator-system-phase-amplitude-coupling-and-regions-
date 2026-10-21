function lam2 = tg_floquet2(g, R, p)
%TG_FLOQUET2  Second-moment Floquet exponent of the homogeneous gamma
%   oscillator when theta runs on a fixed-amplitude cycle
%       xTheta = R cos(wT t),  vTheta = -wT R sin(wT t).
%   Gamma: x'' + (eta0 - g vTheta) x' + wG^2 x = 0 (continuous time).
%   lam2 = (2/T) log max|eig(monodromy)|, T = theta period; this is the
%   growth rate of |u|^2 per unit time. lam2 = -eta0 means no net growth
%   beyond baseline damping. R may be a vector.

R = R(:);
wT = 2*pi*p.fTheta; wG = 2*pi*p.fGamma; eta0 = p.dGamma*wG;
Tp = 2*pi/wT;
cmax = eta0 + g*wT*max(R);
nsub = max(4000, ceil(3*cmax*Tp));          % keeps c*step < ~0.35 for RK4
s = Tp/nsub;
nR = numel(R);
% state: 2x2 fundamental matrix per R, stored as columns [Y11 Y21 Y12 Y22]
Y = repmat([1 0 0 1], nR, 1);
logs = zeros(nR, 1);
t = 0;
for k = 1:nsub
    k1 = rhs(Y, t);
    k2 = rhs(Y + s/2*k1, t + s/2);
    k3 = rhs(Y + s/2*k2, t + s/2);
    k4 = rhs(Y + s*k3, t + s);
    Y = Y + s/6*(k1 + 2*k2 + 2*k3 + k4);
    t = t + s;
    if mod(k, 50) == 0
        n = max(abs(Y), [], 2);
        Y = Y./n; logs = logs + log(n);
    end
end
lam2 = zeros(nR, 1);
for i = 1:nR
    M = [Y(i, 1) Y(i, 3); Y(i, 2) Y(i, 4)];
    lam2(i) = 2*(log(max(abs(eig(M)))) + logs(i))/Tp;
end

    function dY = rhs(Y, t)
        c = eta0 + g*wT*R*sin(wT*t);         % eta0 - g*vTheta
        % scaled coordinates u = (wG x, v): A = [0 wG; -wG -c]
        dY = [wG*Y(:, 2), -wG*Y(:, 1) - c.*Y(:, 2), wG*Y(:, 4), -wG*Y(:, 3) - c.*Y(:, 4)];
    end
end
