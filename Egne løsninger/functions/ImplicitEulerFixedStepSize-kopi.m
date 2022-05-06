function [T,X] = ImplicitEulerFixedStepSize(funJac,ta,tb,N,xa,varargin)
    % Compute step size and allocate memory
    dt = (tb-ta)/N;
    nx = size(xa,1);
    X = zeros(nx,N+1);
    T = zeros(1,N+1);

    tol = 1.0e-8;
    maxit = 100;
    % Eulers Implicit Method
    T(:,1) = ta;
    X(:,1) = xa;
    for k=1:N
        [f, J] = feval(funJac,T(k),X(:,k),varargin{:});
        T(:,k+1) = T(:,k) + dt;
        xinit = X(:,k) + f*dt;
        X(:,k+1) = NewtonsMethodODE(funJac,...
            T(:,k), X(:,k), dt, xinit, tol, maxit, varargin{:});
    end
    % Form a nice table for the result
    T=T';
    X=X';
end