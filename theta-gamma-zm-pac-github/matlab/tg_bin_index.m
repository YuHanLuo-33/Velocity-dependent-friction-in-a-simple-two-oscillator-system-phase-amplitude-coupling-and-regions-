function b = tg_bin_index(phase, nBins)
%TG_BIN_INDEX  Theta-phase bin 1..nBins for phases in [-pi, pi].
if nargin < 2, nBins = 18; end
b = floor((phase + pi)/(2*pi/nBins)) + 1;
b(b > nBins) = nBins;
b(b < 1) = 1;
end
