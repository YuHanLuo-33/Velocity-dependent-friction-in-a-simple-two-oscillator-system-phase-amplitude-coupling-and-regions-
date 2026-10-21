function z = tg_analytic(x)
%TG_ANALYTIC  FFT analytic signal along the first dimension (columns).
n = size(x, 1);
m = zeros(n, 1);
if mod(n, 2) == 0
    m([1, n/2 + 1]) = 1;
    m(2:n/2) = 2;
else
    m(1) = 1;
    m(2:(n + 1)/2) = 2;
end
z = ifft(fft(x).*m);
end
