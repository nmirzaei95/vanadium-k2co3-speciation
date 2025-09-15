function kv = calv_kv(T,I)
% Function for calculating the rate constant for reaction:
% CO2 + HVO4^2- = HVO4CO2^2-.

% Input:    T, temperature (K), scalar
%           Ionic strength (mol/l), scalar

% Output:   kv, rate constant (CO2+HVO4=VC1), (m^3/mol/s), scalar


%% Parameters
E_v = 35.94;                                                % Activation energy (kJ/mol)
lnA_v = 13.24;                                              % Logarithm of Arrhenius-prefactor (Av: m^3/mol/s)
bv_par = [0.23788; -7.1978e-4; -26.146; 7.290e-7];
bv = [1 T 1/T T^2]*bv_par;                                  % Ionic strength pre-factor (m^3/mol)

I = I*1000;                                                 % Conversion: mol/l to mol/m^3

%% Rate constant
ln_kv_inf = lnA_v - E_v*1000/8.314/T;                       % Logarithm of infinite dilution rate constant (kv_inf: m^3/mol/s)
        
ln_kv = ln_kv_inf + bv*(I);                                 % Logarithm of rate constant
kv = exp(ln_kv);                                            % Rate constant (m^3/mol/s)
end