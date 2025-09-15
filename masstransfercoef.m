function [kL] = masstransfercoef(T)
% Function for calculating the liquid-side mass transfer coefficent (kL).
% The kL was experimentally determined in water and its recalculation
% assumes also conditions of water.
% Warning: the liquid-side mass transfer coefficient is setup-dependent.
% However, for reactions studied in this work which are sufficiently fast,
% kL has limited effect on the final results.

% Input:    T, temperature (K), scalar

% Output:   kL, liquid-side mass transfer coefficient (m/s), scalar

%% Constants
d = 0.05;               % Length of stirrer (m)
RPM = 500;              % rotation speed (1/min)
w = RPM/60;             % rotation velocity (1/s)

%% Physical parameters
D = 2.35e-6*exp(-2199./T);                                                          % diffusivity of CO2 in water (m2/s) [1]


A(1) = 1.023096;
A(2) = 0.669483;
A(3) = 640.884;
A(4) = 124.8677;
A(5) = 4694.832*1e-8;

mu = A(5)*exp( A(1)*((A(3)-T)./(T-A(4))).^(1/3) + A(2)*( ((A(3)-T)./(T-A(4))).^(4/3) ) );           % viscosity of water (kg/m/s) [2]


B(1) = -13.418392;
B(2) = 0.6884103;
B(3) = -2.44970115e-3;
B(4) = 3.7060667e-6;
B(5) = -2.11062995e-9;
B(6) = -1.12273895e-13;

rho = 18.015*(B(1) + B(2)*T + B(3)*T.^2 + B(4)*T.^3 + B(5)*T.^4 + B(6)*T.^5);                       % density of water (kg/m3) [2]


% References
% [1] G. F. Versteeg and W. P. M. Van Swaaij, Solubility and diffusivity of
% acid gases (carbon dioxide, nitrous oxide) in aqueous alkanolamine
% solutions, Journal of Chemical & Engineering Data 33, 29 (1988). 

% [2] J. Gmehling, M. Kleiber, B. Kolbe, and J. Rarey, Chemical thermodynamics for process simulation:
% second, completely revised and enlarged edition (John Wiley & Sons, 2019).


%% Dimensionless numbers
Sc = mu./rho./D;                            % Schmidt number
Re = rho.*d^2*w./mu;                        % Reynolds number
Sh = 0.142*Re^0.7*Sc^(1/3);                 % Sherwood number

%% Mass transfer coefficient
kL = Sh*D/d;



end