function [T,X] = ClassicalRungeKuttaFixedStep(...
    fun,tspan,x0,h,varargin)
disp('Lets go')

% Integration interval
t0 = tspan(1);
tf = tspan(2);

% Initial conditions
t = t0;
x = x0;

T = t;
X = x';

%% Algo
while t < tf
    if (t+h > tf)
        h = tf-t;
    end
    % f at x_n
    f = feval(fun, t,x,varargin{:});

    % Take step of size h to obtain x_{n+1}
    [t, x] = ClassicalRungeKuttaStep(...
        fun,t,x,f,h,varargin{:});

     T = [T; t];
     X = [X; x'];
end