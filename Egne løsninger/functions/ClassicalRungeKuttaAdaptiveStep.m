function [T,X,H, fevals] = ClassicalRungeKuttaAdaptiveStep(...
    fun,tspan,x0,h0,abstol,reltol,varargin)
disp('Lets go')
epstol = 0.8;
fevals = 0;
kpow = 0.2; % 1/(order+1)
facmin = 0.1;
facmax = 5;

% Integration interval
t0 = tspan(1);
tf = tspan(2);

% Initial conditions
t = t0;
h = h0;
x = x0;

H = h;
T = t;
X = x';

%% Algo
while t < tf
    if (t+h > tf)
        h = tf-t;
    end
    f = feval(fun, t,x,varargin{:});

    AcceptStep = false;
    while ~AcceptStep
        % Take step of size h
        [t1, x1] = ClassicalRungeKuttaStep(...
            fun,t,x,f,h,varargin{:});

        % Take step of size h/2
        hm=0.5*h;
        [tm, xm] = ClassicalRungeKuttaStep(...
            fun,t,x,f,hm,varargin{:});
        
        fm = feval(fun,tm,xm,varargin{:});
        [t1hat, x1hat] = ClassicalRungeKuttaStep(...
            fun,tm,xm,fm,hm,varargin{:});

        % Error estimation
        e = x1hat-x1;
        r = max(abs(e) ./ max(abstol, abs(x1hat) .*reltol));

        AcceptStep = (r<=1.0);
        if AcceptStep
            t = t+h;
            x = x1hat;
              
            T = [T;t];
            X = [X;x'];
            H = [H;h]; % Save taken step sizes
        end
        % Asymptotic step size controller
        h = max(facmin,min((epstol/r)^kpow,facmax))*h;
        %h = max(facmin, min(sqrt(epstol/r), facmax))*h;
        fevals = fevals + 10; % 10 fevals per full iteration
    end
end