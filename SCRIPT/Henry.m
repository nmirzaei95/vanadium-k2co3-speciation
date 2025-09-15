function H = Henry(T,alfa,CK)
% Function calculating the Henry constant of vanadium promoted PC, based on
% model of Weisenberger-Schumpe. The contribution of vanadate species are
% ignored in the model.

% Inputs:   T, temperature (K), scalar
%           alfa, carbonate conversion (-), scalar
%           CK, potassium concentration (mol/l), scalar

% Outputs:  H, Henry (mol/m3/Pa), scalar


%% Constants


h0_N2O = -8.5e-3;                       % l/mol, [1]
h_CO3 = 0.1423;                         % l/mol, [1-2]
h_HCO3 = 0.0967;                        % l/mol, [1]
h_K = 0.0971;                           % l/mol, [2]
hT_N2O = -0.00001809;                   % l/mol/K, [2]

h_G = h0_N2O + hT_N2O*(T-298.15);       % l/mol

% References:
% [1] S. Weisenberger and A. Schumpe, Estimation of gas solubilities in
% salt solutions at temperatures from 273 K to 363 K, AIChE Journal 42, 298 (1996).

% [2] H. Knuutila, O. Juliussen, and H. F. Svendsen, Density and N2O
% solubility of sodium and potassium carbonate solutions in the temperature
% range 25 to 80 C, Chemical Engineering Science 65, 2177 (2010).


%% Concentration estimations
cCO3 = CK/2*(1-alfa);                           % Approx. concentation of CO3 (mol/l)
cHCO3 = CK*alfa;                                % Approx. concentration of HCO3 (mol/l)



%% Henry constant

H_CO2_w = 3.54 * 1e-7 * exp(2044./T);           % Henry of CO2 in water, (mol/m3/Pa)

% Reference:
% [3] G. F. Versteeg and W. P. M. Van Swaaij, Solubility and diffusivity of
% acid gases (carbon dioxide, nitrous oxide) in aqueous alkanolamine
% solutions, Journal of Chemical & Engineering Data 33, 29 (1988).


S = (h_K + h_G).*(CK) + (h_CO3 + h_G).*cCO3 + (h_HCO3 + h_G).*cHCO3;          % Summation term for the Weisenberger-Schumpe model [1-2]

a = 10.^-S;                     % ratio of Henry of N2O in electrolyte to that of water.
                                % The negative sign arises from the inverse Henry definition of Knuutila

H = H_CO2_w.*a;                 % CO2-N2O analogy. See ref. [3]

end
