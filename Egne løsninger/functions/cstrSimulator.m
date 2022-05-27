function [T,X] = cstrSimulator(method,TimeAx,F,x0,model,abstol,reltol,h0)

T = [];
X = [];

for k = 1:(length(TimeAx)-1)
    tspan = [TimeAx(k), TimeAx(k+1)];
    Fk = F(k);

    %[Tk,Xk,info] = ODEsolver(@cstr3d,tspan,x0,'DOPRI54AdaptiveStep',options,Fk);
    switch upper(model)
        case 'CSTR1D'
            %[Tk,Xk] = ode45(@cstr1d,tspan,x0,options,Fk);
            %[Tk,Xk,info] = ODEsolver(@cstr1d,tspan,x0,'DOPRI54AdaptiveStep',abstol,reltol,Fk);
            solver = ERKSolverErrorEstimationParameters(method);
            [t,x,~] = AdaptiveERKSolverErrorEstimation(@cstr1d,tspan,x0,h0, ...
                solver,abstol,reltol,Fk);

        case 'CSTR3D'
            solver = ERKSolverErrorEstimationParameters(method);
            [t,x,~] = AdaptiveERKSolverErrorEstimation(@cstr3d,tspan,x0,h0, ...
                solver,abstol,reltol,Fk);

        otherwise
            error('Invalid model')
    end
    T = [T;t];
    X = [X;x];

    x0 = x(end,:)';
    %x0 = Xk(:,end);
    fprintf('Iteration %.f of %.f done \n',k,length(F));
end

end