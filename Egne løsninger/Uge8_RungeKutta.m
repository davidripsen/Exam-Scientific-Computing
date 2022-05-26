% Uge8_RungeKutta - Classical Runge Kutta
solver = ERKSolverErrorEstimationParameters('RK44');

% a) Describe Classical Runge-Kutta
% A 4-step method wich is a weighed average of in-between steps when
% estimating x_k+1 from x_k.

%% Test equation
lambda = 0.5;
tspan = [0, 10];
x0 = 1;
h = 1/100;
[Tout,Xout,Eout] = ERKSolverErrorEstimation(@testeq,tspan,x0,h,solver,lambda);

plot(Tout, Xout)
hold on
plot(Tout, exp(lambda*Tout))
legend()
hold off



%% Test equation - Order and Stability
% Order
%  ========================================================================
x0 = 1;
Tf = 50;
Ns = [10 100 1000 10^4 10^5];
hs = (Tf-0)./Ns;
E = zeros(length(hs),1);

for i=1:length(hs)
    [Tout,Xout,Eout] = ...
        ERKSolverErrorEstimation(@testeq,[0 Tf],x0, hs(i),...
                                    solver,lambda)
    E(i) = Xout(end) - x0*exp(Tf*lambda);
end

loglog(hs, abs(E))
xlabel('h')
ylabel('Global error')
polyfit(log(hs), log(abs(E)), 1)
%plot(hspan, E)

%% Van Der Pol
tspan = [0, 50];
mu = 10;
x0 = [0.5; 0.5];
h0 = 1/100; % Initial step size
[Tout,Xout,Eout] = ERKSolverErrorEstimation(@fJacVanDerPol,tspan,x0,h,solver,mu);

plot(Xout(:,1), Xout(:,2))
hold on
[TOUT,YOUT] = ode45(@fJacVanDerPol,tspan,x0,[], mu)
plot(YOUT(:,1), YOUT(:,2))
legend()
hold off
