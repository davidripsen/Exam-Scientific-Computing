function [f, Jac] = fJacVanDerPol(t,x,mu)
    if nargout >= 2 % n arg out
        Jac = zeros(2,2);
        Jac(2,1) = -2*mu*x(1)*x(2)-1.0;
        Jac(1,2) = 1.0;
        Jac(2,2) = mu*(1-x(1)*x(1));
    end
    
    f = zeros(2,1);
    f(1,1) = x(2);
    f(2,1) = mu*(1-x(1)^2)*x(2) - x(1);
     %mu = 10;
     %x0 = [2.0; 0.0];
end