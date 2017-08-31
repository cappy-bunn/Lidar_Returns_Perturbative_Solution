function t = AtmosTransmission(extm,exta)
% Where extm (row vector) is molecular extinction
% and exta (row vector) is aerosol extinction

ncol = size(extm,2);                % Number of columns in either extinction vector
for i=1:ncol
    exti = extm(:,1:i)+exta(:,1:i); % Total extinction
    int_ext(:,i) = trapz(exti);     % Integrate total extinction from 0 to r (keep range info)
end

t = exp(-int_ext);                  % Atmospheric Transmission
end