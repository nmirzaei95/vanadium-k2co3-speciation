function [D_CO2] = diffusivity(T,m)
% Function for calculating the diffusivity of CO2 in aqueous K2CO3. The
% Stokes-Einsten equation modified by Joosten and Danckwerts is used for
% this.

% Input:    T, temperature (K), scalar
%           m, molality (mol K2CO3/kg H2O), scalar


%% Liquid viscosity [1]
f = [2.4251e-1, 4.1839e-2, 6.1826e-3, 1.3116e-3;
    4.7185e-3, -1.1967e-3, 3.3011e-4, -7.6149e-5;
    -3.6792e-5, -1.8959e-6, -1.7439e-6, 8.3178e-7;
    1.6353e-8, 1.8665e-7, -3.1156e-8, -1.4828e-9];

f_i = [0,1,2,3].*ones(4,4);     f_j = [0;1;2;3].*ones(4,4);

mu_rel = 1 + sum(sum(f.*m.^(f_i+1).*(T-273.15).^f_j));
    % ratio of viscosity of electrolyte to viscosity in water


% Reference:
% [1] F. A. Gon¸calves and J. Kestin, The viscosity of Na2CO3 and K2CO3
% aqueous solutions in the range 20-60 C, International Journal of
% Thermophysics 2, 315 (1981).


%% Diffusivity
D_w = 2.35e-6*exp(-2119/T);               % diffusivity of CO2 in water (m2/s) [2]
D_CO2 = D_w*(1/mu_rel)^0.818;               % diffusvity of CO2 in electrolyte (m2/s) [3]


% References:
% [2] G. F. Versteeg and W. P. M. Van Swaaij, Solubility and diffusivity of
% acid gases (carbon dioxide, nitrous oxide) in aqueous alkanolamine
% solutions, Journal of Chemical & Engineering Data 33, 29 (1988). 

% [3] G. E. H. Joosten and P. V. Danckwerts, Solubility and diffusivity of
% nitrous oxide in equimolar potassium carbonate-potassium bicarbonate
% solutions at 25 C and 1 atm, Journal of Chemical and Engineering Data 17,
% 452 (1972).

end