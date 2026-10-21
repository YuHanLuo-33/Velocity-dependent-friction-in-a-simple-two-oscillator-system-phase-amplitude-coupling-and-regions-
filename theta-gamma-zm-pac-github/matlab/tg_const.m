function c = tg_const()
%TG_CONST  Constants of the paper's protocol.
c.fTheta = 4;  c.fGamma = 60;
c.wT = 2*pi*c.fTheta;  c.wG = 2*pi*c.fGamma;
c.h = 5e-4;                         % time step (s)
c.sigma = 200*sqrt(c.h);            % noise intensity, 4.4721
c.nBins = 18;                       % phase bins of the Tort index
end
