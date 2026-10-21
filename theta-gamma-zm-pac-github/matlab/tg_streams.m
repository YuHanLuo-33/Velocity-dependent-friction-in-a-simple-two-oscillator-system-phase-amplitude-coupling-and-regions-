function [z0, streams, inv] = tg_streams(seeds)
%TG_STREAMS  Random-number streams shared by all simulators.
%   For each distinct seed s: initial state from RandStream('mt19937ar', s + 1000003)
%   (randn(4,1) -> [x_theta x_gamma v_theta v_gamma] in units of the uncoupled
%   stationary standard deviations) and step noise from RandStream('mt19937ar', s)
%   (randn(2, n): row 1 theta, row 2 gamma), so a seed gives the same theta path
%   in every simulator.
%   z0 is 4 x numel(seeds); records that share a seed share initial state and noise.
seeds = seeds(:).';
[us, ~, inv] = unique(seeds);
inv = inv(:).';
nU = numel(us);
z = zeros(4, nU);
streams = cell(1, nU);
for k = 1:nU
    s0 = RandStream('mt19937ar', 'Seed', us(k) + 1000003);
    z(:, k) = randn(s0, 4, 1);
    streams{k} = RandStream('mt19937ar', 'Seed', us(k));
end
z0 = z(:, inv);
end
