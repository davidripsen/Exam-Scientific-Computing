function [x, Flag_Diverge] = InexactNewtonsMethodODE(FunJac, tk, xk, dt, xinit, tol, maxit, varargin)
    k=0;
    t=tk+dt;
    x = xinit;
    [f,J] = feval(FunJac,t,x,varargin{:});
    R=x-f*dt-xk;
    I = eye(length(xk));
    Flag_Diverge = false;
    %alphaarray = [];
    while( (k < maxit) && (norm(R,'inf') > tol) )
        normRxk = norm(R);

        k=k+1;
        dRdx=I-J*dt;
        dx = dRdx\R;
        x=x-dx;
        f = feval(FunJac,t,x,varargin{:});
        R=x-dt*f-xk;

        normRxk1 = norm(R);
        alpha = normRxk1 / normRxk;
        if alpha > 1
            Flag_Diverge = true;
            return;
        end
        %alphaarray = [alphaarray, alpha];
    end
end