function f = cstr3d(t,X,F)

% Unpack states 
CA = X(1,:);
CB = X(2,:);
T = X(3,:);

% Parameters
k0 = exp(24.6); % L/(mol*s)
EaDivR = 8500;  % K
CAin = 1.6/2;   % mol/L
CBin = 2.4/2;   % mol/L
Tin = 273.65;   % K
V = 0.105;      % L
dHr = -560;     % kJ/mol
rho = 1.0;      % kg/L
cP = 4.186;     % kJ/(kg*K)

% Rate constant
kT = k0*exp(-EaDivR/T);  

% Reaction rates
r = kT*CA*CB;

beta = -dHr/(rho*cP);

RA = -r;
RB = -2*r;
RT = beta*r;

% Compute derivatives
dCA = F/V*(CAin - CA) + RA;
dCB = F/V*(CBin - CB) + RB;
dT = F/V*(Tin - T) + RT;

% Prepare output
f = [dCA;dCB;dT];

end

