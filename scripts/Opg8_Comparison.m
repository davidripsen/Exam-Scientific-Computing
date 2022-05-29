%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 8 - Discussion and Conlusion
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';
fig.Position(3:4)=[550,400]*1.7;

%% 7.5 - for test equation

lambda = -1;
x0 = 1.0;

abstol = 1e-04;
reltol = abstol;
figure('Position', [100, 100, 1200, 600]);
tspan = [0, 5];
h0 = 1/100; % Initial step size

% Explicit Euler
[T, X] = ExplicitEulerAdaptiveStep(...
    @testeq, tspan,x0,h0,abstol,reltol,lambda);
plot(T, X, DisplayName=sprintf('Explicit Euler,  tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, length(T)-1)) ; hold on

% Implicit Euler
[T,X,iter] = ImplicitEulerAdaptiveStep(...
    @testeq,tspan,x0,h0,abstol,reltol,lambda);
plot(T, X, DisplayName=sprintf('Implicit Euler,  tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on

% Classical Runge-Kutta
[T, X, H, fevals] = ClassicalRungeKuttaAdaptiveStep(...
    @testeq,tspan,x0,h0,abstol,reltol,lambda);
plot(T, X, DisplayName=sprintf('Classical RK,  tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, fevals)) ; hold on


% DOPRI54
solver = ERKSolverErrorEstimationParameters('DOPRI54');
[T,X,E] = AdaptiveERKSolverErrorEstimation(@testeq,tspan,x0,h0, ...
        solver,abstol,reltol,lambda);
plot(T, X, DisplayName=sprintf('DOPRI54, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on

% ESDIRK23
Method = "myESDIRK23";
%[T, X, Gout,info,stats] = ESDIRK(@testeq,@testeqJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,lambda);
%plot(T, X, DisplayName=sprintf('ESDIRK23, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, info.nFun)) ; hold on
%shg % Show graph window

%%% 4b compare with ode45 and ode15
options = odeset('RelTol',reltol,'AbsTol',abstol);
sol45=ode45(@testeq,tspan,x0,options,lambda);
T45 = sol45.x;
X45 = sol45.y;
fevals45 = sol45.stats.nfevals;
sol15=ode15s(@testeq,tspan,x0,options,lambda);
T15 = sol15.x;
X15 = sol15.y;
fevals15 = sol15.stats.nfevals;
plot(T45, X45, DisplayName=sprintf('ode45, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T45), fevals45))
plot(T15, X15, DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T15), fevals15))
title(sprintf("Test equation (lambda = %i)", lambda))
xlabel('T')
ylabel('X')
ylim([0,1])
plot(linspace(0,5,1000), exp(lambda*linspace(0,5,1000)), 'DisplayName','Analytical solution')
legend('location', 'southoutside'); hold off

set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '8a.pdf'))
