% Test equation
function [f  df] = testeq(t, x, lambda)
    x0 = 1; % Change
    f = lambda*x;
    
    df = lambda;
    % Analytical solution
    % x_anal = x0 * exp(lambda*t);
end