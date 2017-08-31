% Doppler Broadened Backscatter Lineshape
function d = DopplerBroadBeta(wn0,wn,m,T)
    c = 3.00e8;          % Speed of light [m/s]
    kb = 1.38065e-23;    % Boltzmann constant [J/K]
    d = (((m*c^2)./(8*pi*wn0^2*kb.*T)).^0.5).*exp(-((m*c^2)./(8*pi*wn0^2*kb.*T)).*(wn'-wn0).^2);
end
