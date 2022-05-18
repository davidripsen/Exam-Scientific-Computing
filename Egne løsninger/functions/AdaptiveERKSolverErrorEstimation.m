function [Tout,Xout,Eout] = ...
        AdaptiveERKSolverErrorEstimation(fun,tspan,x0,h,solver, ...
        abstol, reltol, varargin)
% ERKSOLVERERRORESTIMATION  Fixed step size ERK solver with error est.
%s
%                           Solves ODE systems in the form dx/dt = f(t,x)
%                           with x(t0) = x0. 
%
% Syntax:
% [Tout,Xout,Eout]=ERKSolverErrorEstimation(fun,tspan,x0,h,solver,varargin)

% Adaptive stepsize pars
%abstol=1e-05;
%reltol=1e-05;
epstol=0.8;

% Solver Parameters
s  = solver.stages;     % Number of stages in ERK method
AT = solver.AT;         % Transpose of A-matrix in Butcher tableau
b  = solver.b;          % b-vector in Butcher tableau
c  = solver.c;          % c-vector in Butcher tableau
d  = solver.d;
kpow  = solver.o+1;     % kpow = order + 1

h = h; % INIT

% Size parameters
x  = x0;                % Initial state
t  = tspan(1);          % Initial time
tf = tspan(end);        % Final time
N = (tf-t)/h;           % Number of steps
nx = length(x0);        % System size (dim(x))

% Allocate memory
T  = zeros(1,s);        % Stage T
X  = zeros(nx,s);       % Stage X
F  = zeros(nx,s);       % Stage F

% More init
Tout = t;        
Xout = x';
Eout = 0;

% Algorithm starts here
while t < tf
    if (t+h > tf) % Don't go to far
        h = tf-t;
    end

    AcceptStep = false;
    while ~AcceptStep
        %%% Runge-Kutta
        % Stage 1
        T(1)   = t;
        X(:,1) = x;
        F(:,1) = fun(T(1),X(:,1),varargin{:});
        
        % Stage 2,3,...,s
        T(2:s) = t + h*c(2:s);
        for i=2:s
            X(:,i) = x + F(:,1:i-1)*h*AT(1:i-1,i);
            F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
        end

        % Estimate of apprixmation error
        e = F*h*d;

        % Controlling step size h relative to error
        r = max(abs(e) ./ max(abstol, abs(x + F*h*b) .*reltol));
        h = (epstol/r)^(1/kpow) * h;

        AcceptStep = (r<=1.0);
        if AcceptStep
            t = t+h;
            x = x + F*h*b;
        end
    end
    % Save output when AcceptStep == true
    Tout = [Tout; t];
    Xout = [Xout; x'];
    Eout = [Eout; e];
end