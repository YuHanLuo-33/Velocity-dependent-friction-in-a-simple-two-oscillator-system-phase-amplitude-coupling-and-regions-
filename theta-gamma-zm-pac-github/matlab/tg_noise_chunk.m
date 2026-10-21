function nz = tg_noise_chunk(streams, inv, nC)
%TG_NOISE_CHUNK  Next nC step-noise pairs for every record: 2 x R x nC.
nU = numel(streams);
nz = zeros(2, nU, nC);
for k = 1:nU
    nz(:, k, :) = reshape(randn(streams{k}, 2, nC), 2, 1, nC);
end
nz = nz(:, inv, :);
end
