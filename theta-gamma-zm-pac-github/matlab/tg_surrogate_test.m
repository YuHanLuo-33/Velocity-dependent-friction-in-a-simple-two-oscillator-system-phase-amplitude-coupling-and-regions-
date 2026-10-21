function S = tg_surrogate_test(P, nSur, stream)
%TG_SURROGATE_TEST  Within-cycle circular-shift surrogates.
%   Each theta cycle (from phase wraps) of the gamma envelope is shifted
%   circularly by an independent random offset; theta phase is kept. The same
%   offsets are applied to the raw and the detrended envelope. Returns [] if the
%   record has fewer than 8 cycles.
K = 18;
phase = P.phase(:);
N = numel(phase);
wrapped = angle(exp(1i*phase));
if median(diff(unwrap(wrapped))) > 0
    starts = find(diff(wrapped) < -pi) + 1;
else
    starts = find(diff(wrapped) > pi) + 1;
end
a = [1; starts(:)];
b = [starts(:) - 1; N];
L = b - a + 1;
keep = L >= max(4, 0.5*median(L));
a = a(keep); L = L(keep);
if numel(a) < 8
    S = [];
    return
end
onehot = sparse((1:N).', P.bins(:), 1, N, K);
counts = P.counts(:);
rawNull = zeros(1, nSur);
detNull = zeros(1, nSur);
block = 200;
done = 0;
while done < nSur
    nb = min(block, nSur - done);
    idx = repmat((1:N).', 1, nb);
    for j = 1:numel(a)
        Lj = L(j);
        off = randi(stream, Lj, 1, nb) - 1;
        p = (0:Lj-1).';
        idx(a(j):a(j)+Lj-1, :) = a(j) + mod(p - off, Lj);
    end
    rawNull(done+1:done+nb) = tg_mi_from_means(full(onehot.'*P.amp(idx))./counts);
    detNull(done+1:done+nb) = tg_mi_from_means(full(onehot.'*P.det(idx))./counts);
    done = done + nb;
end
S.raw = summ(P.rawMI, rawNull);
S.det = summ(P.detMI, detNull);
S.nCycles = numel(a);
S.nValid = nSur;
end

function s = summ(obs, null)
s.p = (1 + sum(null >= obs))/(numel(null) + 1);
s.delta = obs - median(null);
s.z = (obs - mean(null))/max(std(null), 1e-300);
s.nullMedian = median(null);
end
