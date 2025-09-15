function [CK] = X2CK2(X,T)
% Function for calculating molarity of K2CO3 in the solvent from its mass
% fraction and temperature via its density. Influence of V2O5 on the
% density is neglected.

% Input:    X, mass fraction of K2CO3 equivalent (unpromoted and tht=0), scalar
%           T, temperature (K)

% Output:   CK, potassium concentration (mol/l), scalar

%% Conversion of mass fraction to molarity
Mw_K2CO3 = 138.2;                              % Molecular weight of K2CO3, g/mol
C_m_K2CO3 = X/(Mw_K2CO3*1e-3)./(1-X);          % Solvent molality (mol K2CO3/kg H2O)

%% Liquid density
d = [1.0006e3, 1.2479e2, -1.0544e1, 6.8296e-1, -2.3941e-2;
    -1.4378e-3, -2.3933e-1, 2.0655e-2, 1.911e-3, -2.2461e-4;
    -5.7652e-3, 7.1823e-4, 1.0597e-3, -2.7032e-4, 1.7181e-5;
    1.6986e-5, 1.0135e-5, -1.2824e-5, 2.7504e-6, -1.6771e-7];

d_i = [0, 1, 2, 3, 4].*ones(4,5);       d_j = [0; 1; 2; 3].*ones(4,5);
rho_l = zeros(size(X));

for i = 1:size(X,1)
    rho_l(i) = sum(sum(d.*C_m_K2CO3(i).^d_i.*(T(i)-273.15).^d_j));            % Density of liquid, [kg/m^3]
end

% Reference
% F. A. Gon¸calves and J. Kestin, The viscosity of Na2CO3 and K2CO3 aqueous solutions in the range
% 20-60 C, International Journal of Thermophysics 2, 315 (1981).    


Mw_K2CO3 = 138.2;                       % Molecular weight of K2CO3, [g/mol]
C_K2CO3 = X.*rho_l/138.2;               % Molarity of K2CO3 [mol/l]

CK = C_K2CO3*2;                 

end
