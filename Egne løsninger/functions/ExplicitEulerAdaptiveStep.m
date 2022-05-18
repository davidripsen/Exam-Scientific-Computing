function [T,X,fcount,nreject] = ExplicitEulerAdaptiveStep(...
    fun,tspan,x0,h0,abstol,reltol,varargin)
epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

t0 = tspan(1);
tf = tspan(2);
fcount = 0;
nreject = 0;

% Initial condtions
t = t0;
h = h0;
x = x0;

% Output
T = t;
X = x';
%% Main algorithm
while t < tf
    if (t + h >tf)
        h = tf-t;
    end
    fcount = fcount+1;
    f = feval(fun,t,x,varargin{:});

    AcceptStep = false;
    while ~AcceptStep
        %Take step of size h
        x1 = x+h*f;

        %Take step of size h/2
        hm = 0.5*h;
        tm = t + hm;
        xm = x + hm*f;
        fcount = fcount +1;
        fm = feval(fun,tm,xm,varargin{:});
        x1hat = xm + hm*fm;

        % Error estimation
        e = x1hat-x1; % Estimate of global error
        r = max(abs(e) ./ max(abstol,abs(x1hat).*reltol));

        AcceptStep = (r <= 1.0);
        if AcceptStep
            t = t+h;
            x = x1hat;

            T = [T;t];
            X = [X;x'];
        else
            nreject = nreject+1;
        end
        % Asymptotic step size controller
        h = max(facmin, min(sqrt(epstol/r), facmax))*h;
    end
end