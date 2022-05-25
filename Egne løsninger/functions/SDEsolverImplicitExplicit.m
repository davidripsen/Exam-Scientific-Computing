function X=SDEsolverImplicitExplicit(ffun,gfun,T,x0,W,varargin)
tol = 1.0e-8;
maxit = 100;

N = length(T);
nx = length(x0);
X = zeros(nx,N);

X(:,1) = x0;
k=1;
[f,~] = feval(ffun,T(k),X(:,k),varargin{:});
for k=1:N-1
    g = feval(gfun,T(k),X(:,k),varargin{:});
    dt = T(k+1)-T(k);
    dW = W(:,k+1)-W(:,k);
    psi = X(:,k) + g*dW;
    xinit = psi + f*dt;
    [X(:,k+1),f,~] = SDENewtonSolver(...
            ffun,...
            T(:,k+1),dt,psi,xinit,...
            tol,maxit,varargin{:});
end