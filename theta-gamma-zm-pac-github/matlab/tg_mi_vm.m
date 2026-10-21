function m = tg_mi_vm(k, K)
%TG_MI_VM  Continuum Tort index of a von Mises envelope exp(k cos phi):
%   MI(k) = [k I1(k)/I0(k) - ln I0(k)] / ln K   (closed form of the paper).
if nargin < 2, K = 18; end
i0 = besseli(0, k, 1);                 % scaled: I0(k) exp(-k)
i1 = besseli(1, k, 1);
m = (k.*i1./i0 - (log(i0) + k))/log(K);
end
