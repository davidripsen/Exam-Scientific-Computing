% Uge 3 - Euler Solvers + PDEs
% Start med exercises fra Lect2C
t = 1;
x = [2.0; 2.0];
mu = 3;
[f, Jac] = fJacVanDerPol(t, x, mu)

%% EXPLICIT-EULER (Jac Van Der Pol)
x0 = [2; 0];
mu = 3;
    % Solver specifics
N = 10000;
t0 = 0;
tN = 100;

[T, Xex] = ExplicitEulerFixedStepSize(@fJacVanDerPol,t0,tN, N, x0, mu)
plot(Xex(:,1), Xex(:,2))

%% IMPLICIT-EULER (Jac Van Der Pol)
xa = [2; 0];
ta = 0;
tb = 100;
N = 10000;
mu = 3;

[T,Xim] = ImplicitEulerFixedStepSize(@fJacVanDerPol,ta,tb,N,xa,mu)
plot(Xim(:,1), Xim(:,2))

%% Comparisons of solutions
% plot differences
plot(Xex(:,1), Xex(:,2))
hold on
plot(Xim(:,1), Xim(:,2))
hold off

%% %%%%%%%%%%% Plug Flow Reactor (1st order kinetics) %%%%%%%%%%%%%%%%%%%%
% Goal: Model the concentration of A, determined by location in the pipe
% (z) and time (t). z = L is the end of the pipe.
% C = concentration of A
% N = mol of A
% R = production rate of A (negative in this case, because A --> P)
% v =  
dz = 0.1;
k = 0.5;
Nz = 10;
v = 10;
k = 0.8;
C0 = 3;
Cj1 = 4;

pipeflow(0,2,C0,Cj1,Nz,k)

% MANGLER:
% ExplicitEulerFixedStepsize(pipeflow,..., <Cj1 = NewtonsMethod()>)

%% Functions
% Something in the lines of
function Cdot = pipeflow(t,j,C,Cj1,Nz,k)
    % Reaction and production rates
    r = zeros(Nz); r(j) = k*C(j);
    R = zeros(Nz); R(j) = -r(j);
    
    % C(j+1)=
    Cj1 = Cj1; % Fundet ved NewtonsMethodODE() (MANGLER)
    
    % Fluxes
    N = v*C(j) - (0 < j < Nz) * D * (Cj1 - C)/dz;

    % Finally: Put things together in the (Spatial) discretization of PDE
    Cdot = -N + R(j)
        % Noget med noget second order finite difference?
end