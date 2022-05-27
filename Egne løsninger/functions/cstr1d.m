function f = cstr1d(t,X,F)

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

% Unpack state
T = X;

% Rate constant
kT = k0*exp(-EaDivR/T);  

beta = -dHr/(rho*cP);

% Approximation of concentrations of A and B
CA = CAin + 1/beta*(Tin-T);
CB = CBin + 2/beta*(Tin-T);

% Reaction rate
r = kT*CA*CB;

RT = beta*r;

dT = F/V*(Tin-T) + RT;

f = dT;

end