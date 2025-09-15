% The script is written as an addendum to the following studies:
% (1) Experiments and kinetic modeling of absorption rates of CO2 into
% unpromoted K2CO3 solutions at low to high solvent loading
% (2) Kinetic and mechanistic study of CO2 absorption into vanadium-promoted
% aqueous K2CO3, by Mirzaei N. and Babler M. U.

% The collective of the provided script and functions use the models
% introduced in this article to recalculate (among others), species
% concentrations, reaction rate constants, and mass transfer coefficients
% for aqueous K2CO3 with or without vanadium pentoxide.

% For the models and equations used, the user is referred to the
% accompanied documentation and the original articles themselves.

% N. Mirzaei Sep. 2025
% v1


clc;
close all;


%% Inputs
CK = 4.4;                       % concentration of potassium, mol/l (CK = 2[K2CO3])
CV = 0.4;                       % concentration of vanadium, mol/l (CV = 2[V2O5])
tht = 0.2;                      % solvent loading -
T = 333;                        % temperature, K



%% Auxiliary parameters
Mw_K2CO3 = 138.2;               % molecular weight of K2CO3, g/mol
alfa = tht + 2*CV./CK;          % carbonate conversion


%% Physico-chemical parameters


%%%% Equivalent molality and mass fraction of K2CO3 %%%%
X_DB = linspace(0.01,0.4);                      % Database: mass fraction of K2CO3
CK_DB = zeros(size(X_DB));
for i = 1:length(X_DB)
    [CK_DB(i)] = X2CK2(X_DB(i),298);
end

X = interp1(CK_DB,X_DB,CK);                     % Mass fraction of K2CO3
m = X/(Mw_K2CO3*1e-3)./(1-X);                   % Solvent molality (mol K2CO3/kg H2O)


kL = masstransfercoef(T);       % liquid-side mass transfer coefficient (m/s)

%%%% Henry constant %%%%
H = Henry(T,alfa,CK);           % Henry constant (mol/m^3/Pa)

%%%% Diffusivity %%%%
D = diffusivity(T,2.4);         % diffusivity (m^2/s)



%% Species concentrations
lim_DB = [8.5 12];                                                  % pH limits for the database
N = 300;                                                            % database resolution; higher resolution = better results/slower convergence
    % lim_DB and N can be adjusted based on user's needs

pH_DB = linspace(lim_DB(1),lim_DB(2),N);                            % pH database,
gs = CV*0.1*ones(size(pH_DB));                                      % guess for cH2VO4 (0<gs<CV generally works fine)
options = optimoptions('fsolve','StepTolerance',1e-10);
x = fsolve(@ (x) ChEq_V(x,pH_DB,CK,CV,tht,T), gs,options);          % concentration of cH2VO4, mol/l

[~,c,I] = ChEq_V(x,pH_DB,CK,CV,tht,T);
    % c: species concentrations corresponding to CK, CV, tht, and T (mol/l)
    % order [CO3, HCO3, CO2, H, OH, HVO4, H2VO4, VO4, V2O7, HV2O7, H2V2O7, HV3O10 V4O13, V4O12, V5O15, VC1, VC2]
    % I: ionic strength (mol/l)

%% Rate constants

k2 = calc_k2(T,alfa,CK);                        % forward rate constant for reaction CO2 + OH = HCO3 (m^3/mol/s)

kv = calv_kv(T,I);                              % forward rate constant for reaction CO2 + HVO4^2- = HVO4CO2^2- (m^3/mol/s)

cOH = c(5)*1000;                                % concentration of OH^- (mol/m^3)
cHVO4 = c(6)*1000;                              % concentration of HVO4^2- (mol/m^3)

k1 = k2*cOH + kv*cHVO4;                         % pseudo-first order rate constant (1/s)


%% Overall mass transfer coefficient
M = k1.*D/kL.^2;
E = sqrt(M)./tanh(sqrt(M));
Kg = kL.*H.*E;

