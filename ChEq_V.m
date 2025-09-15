function [F,c,I] = ChEq_V(x,pH,CK,CV,tht,T)
% Objective function for estimating the speciation of solvent for given
% macroscopic composition, i.e., potassium, vanadium, and carbon amounts

% Input:    x, concentration of H2VO4- corresponding to pH (mol/l), 1-by-m vector
%           pH, from database (-), 1-by-m vector
%           CV, concentration of vanadium (mol/l), scalar
%           CK, concentration of potassium (mol/l), scalar
%           tht, solvent loading (-), scalar
%           T, temperature (K), scalar

% Output:  F, objective function corresponding to x and pH, 1-by-m vector
%          c, species concentrations, corresponding to CK, CV, tht, T (mol/l), 17-by-1 vector
%               [CO3, HCO3, CO2, H, OH, HVO4, H2VO4, VO4, V2O7, HV2O7, H2V2O7, HV3O10 V4O13, V4O12, V5O15, VC1, VC2]
%          I, ionic strength (mol/l), scalar

% N. Mirzaei Sep. 2025
% v1


cH2VO4 = x;                 % H2VO4 concentration




%% Equilibrium constants                                                                             References                  Stoichiometry
K1 = exp(-1203.01 + 68359.6./T + 188.444*log(T) -0.206424*T - 4712910./T.^2);       % [mol/l]       [1]                     CO2 + H2O = HCO3 + H
K2 = exp(175.36 - 7230.6./T - 30.6509*log(T) + 0.01315*T - 372805./T.^2);           % [mol/l]       [1]                     HCO3 = CO3 + H
Kw = exp(140.932 - 13445.9/T - 22.4773*log(T));                                     % [mol2/l2]     [1]                     H2O = OH + H
Kv1 = 10^-2;                                                                        % [mol/l]       Set to be very low      H3VO4 = H2VO4 + H
Kv2 = exp(-14.999-1738.06/T);                                                       % [mol/l]       [1]                     H2VO4 = HVO4 + H
Kv3 = 10^-21.31/Kv2;                                                                % [mol/l]       [2]                     HVO4 = VO4 + H
Kvc1 = exp(13.342+3863.57/T);                                                       % [l2/mol2]     [1]                     CO3+H2VO4+H = HVO4CO2+2H2O
Kvc2 = exp(26.474+6452.61/T);                                                       % [l4/mol4]     [1]                     2CO3+H2VO4+2H = VO4(CO2)2+2H2O
K3 = exp(-47.72+1083.93/T);                                                         % [mol/l]       [1]                     2H2VO4 = V2O7 + 2H + H2O
K4 = exp(-18.13+251.94/T);                                                          % [-]           [1]                     2H2VO4 = HV2O7 + H + H2O
K5 = exp(0.8778 + 1200.18/T);                                                       % [l/mol]       [1]                     2H2VO4 = H2V2O7 + H2O
K6 = exp(-21.485+2193.35/T);                                                        % [l/mol]       [1]                     3H2VO4 = HV3O10 + H + 3H2O
K7 = exp(-42.167+1484.18/T);                                                        % [l/mol]       [1]                     4H2VO4 = V4O13 + 2H + 3H2O
K8 = exp(-12.058+8369/T);                                                           % [mol3/l3]     [1]                     4H2VO4 = V4O12 + 4H2O
K9 = exp(-25.685+12783/T);                                                          % [l4/moll4]    [1]                     5H2VO4 = V5O15 + 5H2O

% References:
% [1] M. Imle, J. Kumelan, D. Speyer, N. McCann, G. Maurer, and H. Hasse, Solubility of carbon dioxide in
% activated potash solutions in the low and high gas loading regions, Industrial & Engineering Chemistry
% Research 52, 13477 (2013).
% [2] D. C. Crans and A. S. Tracey, The chemistry of vanadium in aqueous and nonaqueous solution, in
% Vanadium Compounds: Chemistry, Biochemistry, and Therapeutic Applications (ACS Publications,
% 1998).



%% Concentrations of H and OH and auxiliary parameters (see documentation)
cH = 10.^-pH;              cOH = Kw./cH;

a = [1 + 2*Kv2./cH + 3*Kv2*Kv3./cH.^2;
    4*K3./cH.^2 + 3*K4./cH + 2*K5;
    4*K6./cH;  6*K7./cH.^2 + 4*K8;    5*K9*ones(size(pH))];


b = [cH/Kv1 + 1 + Kv2./cH + Kv2*Kv3./cH.^2;
    2*K3./cH.^2 + 2*K4./cH + 2*K5;
    3*K6./cH;  4*K7./cH.^2 + 4*K8;    5*K9*ones(size(pH))];

g1 = CK;
g2 = CV;

e1 = 2 + cH/K2;
e3 = 1 + cH/K2 + cH.^2/K1/K2;

l = cH - cOH;




%% concentration of carbonates as function of H2VO4
cCO3 = (g1 - 3*g2 + l - sum( (a-3*b) .* ( cH2VO4.^[1;2;3;4;5] ) ) )./ (e1 - Kvc1*cH.*cH2VO4);


%% Objective function: vanadium balance
lft = g2;
rgt = sum( b .* cH2VO4.^[1;2;3;4;5] ) + Kvc1*cH.*cH2VO4.*cCO3 + Kvc2*cH.^2.*cH2VO4.*cCO3.^2;

F = lft-rgt;





if nargout>1

%% Create database (struct D) 
    D.cH = cH;                                  % concentration of H (mol/l), 1-by-m vector
    D.cOH = cOH;                                % concentration of OH (mol/l), 1-by-m vector
    D.cH2VO4 = cH2VO4;                          % concentration of H2VO4 (mol/l), 1-by-m vector
    D.cCO3 = cCO3;                              % concentration of CO3 (mol/l), 1-by-m vector
    D.cHCO3 = cH/K2.*cCO3;                      % concentration of HCO3 (mol/l), 1-by-m vector
    D.cHVO4 = Kv2./cH.*cH2VO4;                  % concentration of HVO4 (mol/l), 1-by-m vector
    D.cVO4 = Kv2*Kv3./cH.^2.*cH2VO4;            % concentration of VO4 (mol/l), 1-by-m vector
    D.cV2O7 = K3./cH.^2.*cH2VO4.^2;             % concentration of V2O7 (mol/l), 1-by-m vector
    D.cHV2O7 = K4./cH.*cH2VO4.^2;               % concentration of HV2O7 (mol/l), 1-by-m vector
    D.cH2V2O7 = K5*cH2VO4.^2;                   % concentration of H2V2O7 (mol/l), 1-by-m vector
    D.cHV3O10 = K6./cH.*cH2VO4.^3;              % concentration of HV3O10 (mol/l), 1-by-m vector
    D.cV4O13 = K7./cH.^2.*cH2VO4.^4;            % concentration of V4O13 (mol/l), 1-by-m vector
    D.cV4O12 = K8*cH2VO4.^4;                    % concentration of V4O12 (mol/l), 1-by-m vector
    D.cV5O15 = K9*cH2VO4.^5;                    % concentration of V5O15 (mol/l), 1-by-m vector
    D.cVC1 = Kvc1*cH2VO4.*cCO3.*cH;             % concentration of VC1 (mol/l), 1-by-m vector
    D.cVC2 = Kvc2*cH2VO4.*cCO3.^2.*cH.^2;       % concentration of VC2 (mol/l), 1-by-m vector
    D.cCO2 = cH/K1.*D.cHCO3;                    % concentration of (free) CO2 (mol/l), 1-by-m vector

    D.Ccarbon = D.cCO3 + D.cHCO3 + D.cCO2 + D.cVC1 + 2*D.cVC2;              % total carbon concentration (mol/l) corresponding to pH within database
    D.tht = D.Ccarbon / (CK/2) - 1;             % solvent loading corresponding to each pH within the database

%% Interpolate from the database species concentrations corresponding to solvent loading tht
    cCO3 = interp1(D.tht,D.cCO3,tht);
    cHCO3 = interp1(D.tht,D.cHCO3,tht);
    cCO2 = interp1(D.tht,D.cCO2,tht);
    cH = interp1(D.tht,D.cH,tht);
    cOH = interp1(D.tht,D.cOH,tht);
    cHVO4 = interp1(D.tht,D.cHVO4,tht);
    cH2VO4 = interp1(D.tht,D.cH2VO4,tht);
    cVO4 = interp1(D.tht,D.cVO4,tht);
    cV2O7 = interp1(D.tht,D.cV2O7,tht);
    cHV2O7 = interp1(D.tht,D.cHV2O7,tht);
    cH2V2O7 = interp1(D.tht,D.cH2V2O7,tht);
    cHV3O10 = interp1(D.tht,D.cHV3O10,tht);
    cV4O13 = interp1(D.tht,D.cV4O13,tht);
    cV4O12 = interp1(D.tht,D.cV4O12,tht);
    cV5O15 = interp1(D.tht,D.cV5O15,tht);
    cVC1 = interp1(D.tht,D.cVC1,tht);
    cVC2 = interp1(D.tht,D.cVC2,tht);
    
    
    c(1,:) = cCO3;
    c(2,:) = cHCO3;
    c(3,:) = cCO2;
    c(4,:) = cH;
    c(5,:) = cOH;
    c(6,:) = cHVO4;
    c(7,:) = cH2VO4;
    c(8,:) = cVO4;
    c(9,:) = cV2O7;
    c(10,:) = cHV2O7;
    c(11,:) = cH2V2O7;
    c(12,:) = cHV3O10;
    c(13,:) = cV4O13;
    c(14,:) = cV4O12;
    c(15,:) = cV5O15;
    c(16,:) = cVC1;
    c(17,:) = cVC2;

%% Ionic strength
    I = 1*CK + 4*cCO3 + 1*cHCO3 + 0*cCO2 + 1*cH + 1*cOH +...
        1*cH2VO4 + 4*cHVO4 + 9*cVO4 + 16*cV2O7 + 9*cHV2O7 + 4*cH2V2O7 + 16*cHV3O10 + 36*cV4O13 + 16*cV4O12 + 25*cV5O15 +...
        4*cVC1 + 9*cVC2;
    I = I/2;

end    
