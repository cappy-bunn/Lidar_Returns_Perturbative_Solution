% Calculates G1 for a modeled system. Modeled system can be found by
% varying the lapse rate (gamma).

% gamma is the adiabatic lapse rate [K/km]


function d = ModeledG1(gamma)

global wl_on wl_off Ts Ps bg bl blw sa sm c n0 kb h m Epulse tpulse...
    rangebin r rkm n fwidth wn0_on f0_on f1_on f2_on f_on f wn1_on...
    wn2_on wn_on_step wn_on wn0_off f0_off f1_off f2_off wn1_off wn2_off...
    wn_off_step wn_off; %#ok<NUSED>

% --------- Profiles ---------
T = Ts+gamma.*rkm;                      % Temperature [K]
P = Ps*(Ts./T).^(-5.2199);              % Pressure [atm]
Pkpa = 101.325*P;                       % Pressure [kPa]
nO2 = 0.2096*n0*T.*P/(273.15*1);        % O2 number density [m^-3]

betam1 = 374280*Pkpa./(T*wl_on^4);      % ONLINE Molecular Backscatter

% Aerosol Backscatter
k = size(r);
k = k(2);
for j=1:1:k
    range = rkm(j);
    betaa(j) = bg/100;
    if range<bl+blw/2
        betaa(j) = bg*(1-(bl-blw/2)/10)*(1-(range-(bl-blw/2))/blw);
    end
    if range<bl-blw/2
        betaa(j) = bg*(1-range/10);
    end
end
exta = sa*betaa;                      % Aerosol Extinction (Casey ~Eq. 45)

% ONLINE Doppler Broadened Backscatter Lineshape
D1_unnorm = DopplerBroadBeta(wn0_on,wn_on,m,T);
D1_norm_func = trapz(D1_unnorm.*wn_on_step);
D1 = D1_unnorm./D1_norm_func./2;                   % Normalized lineshape

% -------- Absorption Cross Section Calculation ---------
[~,sigmaV1] = AbsorptionCrossSection(wl_on,wl_on,T,P,r);

alpham1 = sigmaV1.*repmat(nO2',1,301)/(100^2);  % repeat array nO2 to create matrix same size as sigmaV1
% alpham is absorption coefficient for given height and frequency (within our range)

% -------- Molecular Absorption ---------
% First find optical depth from absorption coefficient (alpham).
optdepth1 = OpticalDepth(alpham1,rangebin);  % Optical depth for online

% Next find molecular absorption: Tm = exp[-int(alpha)dr from r0 to r]
Tm1 = exp(-optdepth1);       % Fraction of online light reaching a range r

% ------ Return Lineshape ------
betat1 = betaa+betam1;                      % Total online backscatter

aero_back1 = repmat((betaa./betat1),301,1); % Aerosol reflects delta function lineshape
aero_back1(1:150,:)=0;
aero_back1(152:301,:)=0;

gi1 = aero_back1+(betam1./betat1).*D1;      % Total return lineshape online

% ---- Etalon Transmission ----
wl_start = wl_on-0.01972;           % Beginning of wavelength range [nm] (equivalent to 10 GHz away from online)
wl_end = wl_on+0.01972;             % End of wavelength range [nm]
wl_step = (wl_end-wl_start)/(n-1);  % Wavelength step size (made so the range has 301 points)
lambda = wl_start:wl_step:wl_end;   % Wavelength range

R = 0.92956;                            % Etalon reflectivity
FSR = wl_off-wl_on;                     % Free Spectral Range (FSR=c/(2*nL)) difference between online and offline wl [nm]
theta = (2*pi*(lambda-wl_on))/FSR;                      % Round trip phase accumulation
etalonT_unnorm = (1-R)^2./(pi*(1+R^2-2*R*cos(theta)));  % Etalon transmission (before normalization)
etalon_norm_factor = etalonT_unnorm(:,151);             % Normalization factor
etalonT = etalonT_unnorm/etalon_norm_factor;            % Normalized etalon transmission

% ----- Calculate G1 ------
% Find the derivative wrt r of gi1:
num_col = size(gi1,2);
ii=1:num_col-2;
jj=3:num_col;
for ix = 1:length(ii)
    ii(ix);
    jj(ix);
    dgi1_dr_new(:,ii) = (gi1(:,jj)-gi1(:,ii))/(2*rangebin);
end

% Repeat the first and last columns to maintain original matrix dimensions:
dgi1_dr = [dgi1_dr_new(:,1) dgi1_dr_new dgi1_dr_new(:,37)];

int_numG1 = trapz(dgi1_dr.*etalonT'.*Tm1');
int_denG1 = trapz(gi1.*etalonT'.*Tm1');
d = int_numG1./int_denG1;         % Correction factor G1


end
