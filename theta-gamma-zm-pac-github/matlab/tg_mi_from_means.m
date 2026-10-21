function mi = tg_mi_from_means(m)
%TG_MI_FROM_MEANS  Tort MI from mean amplitudes per bin (bins along dimension 1).
K = size(m, 1);
p = m./sum(m, 1);
mi = (log(K) + sum(p.*log(p), 1))/log(K);
mi = min(max(mi, 0), 1);
end
