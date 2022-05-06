function [T,X,iter] = ImplicitEulerAdaptiveInexactNewt(...
    funJac,tspan,x0,h0,abstol,reltol,varargin)
epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

t0 = tspan(1);
tf = tspan(2);

iter = 0;
% Initial condtions
t = t0;
h = h0; % = dt (step size - bliver modificeret)
x = x0;

% Output
T = t;
X = x';

% Newton stuff
tol = 1e-05;
maxit = 1000;

%% Main algorithm
while t < tf
    iter = iter +1;
    if (t + h >tf)
        h = tf-t;
    end
    f = feval(funJac,t,x,varargin{:});
    AcceptStep = false;
    while ~AcceptStep
        %Take step of size h
        xinit = x + h*f; % Gæt på x+1
        [x1, Flag_Diverge] = InexactNewtonsMethodODE(funJac, t, x, h, xinit, tol, maxit, varargin{:});
        if ~Flag_Diverge
            %Take step of size h/2
            hm = 0.5*h;
            xinitm = x + hm*f;
            tm = t + hm;
    
            xm = InexactNewtonsMethodODE(funJac, tm, x, hm, xinitm, tol, maxit, varargin{:});
            
            xinitm2 = x + 2*hm*f ; % = xinit (og derfor redundant)
            x1hat = InexactNewtonsMethodODE(funJac, tm, xm, hm, xinitm2, tol, maxit, varargin{:});
    
            % Error estimation
            e = x1hat-x1; % Estimate of global error
            r = max(abs(e) ./ max(abstol,abs(x1hat).*reltol));
    
            AcceptStep = (r <= 1.0);
            if AcceptStep
                t = t+h;
                x = x1hat;
    
                T = [T;t];
                X = [X;x'];
            end
            % Asymptotic step size controller
            h = max(facmin, min(sqrt(epstol/r), facmax))*h;
        else
            h=h/2;
        end
    end
end