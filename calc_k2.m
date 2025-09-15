function k2 = calc_k2(T,alfa,CK)
% Function for calculating the rate constant for reaction:
% CO2 + OH- = HCO3-.
% Warning: Model applies to case of alfa>0!

% Input:    T, temperature (K), scalar
%           alfa, carbonate conversion (-), scalar
%           CK, potassium concentration (mol/l), scalar

% Output:   k2, rate constant (CO2+OH=HCO3), (m^3/mol/s), scalar


%% Equilibrium constants
K2 = exp(175.36 - 7230.6./T - 30.6509*log(T) + 0.01315*T - 372805./T.^2);           % [mol/l]       HCO3 = CO3 + H
Kw = exp(140.932 - 13445.9/T - 22.4773*log(T));                                     % [mol2/l2]     H2O = OH + H

% Reference:
% [1] M. Imle, J. Kumelan, D. Speyer, N. McCann, G. Maurer, and H. Hasse, Solubility of carbon dioxide in
% activated potash solutions in the low and high gas loading regions, Industrial & Engineering Chemistry
% Research 52, 13477 (2013).

%% Ion-contribution model coefficients
A1 = [0.5908e0; -0.1648e-2; 1.5212e-6; -6.9966e1];              % Parameters for bta1
A2 = [1.08455e2; -3.1289e-1; 3.0092e-4; -1.2535e4];             % Parameters for bta2

bta1 = [ones(size(T)) T T.^2 1./T]*A1;            bta2 = [ones(size(T)) T T.^2 1./T]*A2;            % bta parameters from unpromoted solvent

% Reference:
% N. Mirzaei, F.Walthert, E. Kantarelis, and M. U. Babler, Experiments and
% kinetic modeling of absorption rates of CO2 into unpromoted K2CO3
% solutions at low to high solvent loading, Separation and Purification
% Technology, 134622 (2025)


%% Hydroxide concentrations at reference and at alfa
cOH_ref = -1/2*Kw./K2 + 1/2*sqrt( 2*CK.*Kw./K2 + (Kw./K2).^2 );                                 % Hydroxide conc. at reference alfa=0, mol/l
cOH = -1/2*CK.*alfa -1/2*Kw./K2 + 1/2*sqrt( (CK.*alfa).^2 + 2*CK.*Kw./K2 + (Kw./K2).^2 );       % Hydroxide conc. at alfa, mol/l

cOH_ref = cOH_ref*1000;
cOH = cOH*1000;                                                 % Conversion: mol/l to mol/m3


%% Rate constants

% at reference point (alfa=0)

%%%% 25 wt% K2CO3 %%%%
k2_ref = exp(22.01 - 5771./T);              % Ref: [1]



%%%% not 25 wt% K2CO3 %%%%
lnA = 0.24*3/2*CK + 26.4 - log(1000);
Ea = 47.03*1000;
k2_ref = exp(lnA - Ea/8.314/T);             % Ref: [2]

% References:
% [1] N. Mirzaei, F.Walthert, E. Kantarelis, and M. U. Babler, Experiments and
% kinetic modeling of absorption rates of CO2 into unpromoted K2CO3
% solutions at low to high solvent loading, Separation and Purification
% Technology, 134622 (2025)

% [2] X. Ye and Y. Lu, Kinetics of CO2 absorption into uncatalyzed
% potassium carbonate-bicarbonate solutions: Effects of CO2 loading and
% ionic strength in the solutions, Chemical Engineering Science 116, 657
% (2014).


k2 = k2_ref .* exp( bta1.*( (CK*1000) ).*alfa + bta2.*(cOH - cOH_ref) );            % Rate constant at desired point, m^3/mol/s
end