%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 2 - Explicit ODE solver
clc;
plotpath = '../plots/';

%% 1. Describe the Explicit Euler
disp('Bla bla bla')

%% 2. 
% se ExplicitEulerFixedStepSize

%% 3.
% se ExplicitEulerAdaptievStepSize

%% 4. Adaptive on Jac Van Der Pol
tspan = [0, 40];
mu = 3; % mu = 20
x0 = [1.0; 1.0];
h0 = 1/100; % Initial step size
abstol = 1e-06;
reltol = 1e-06;
[T, X] = ExplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

%%% 4b compare with ode45 and ode15
options = odeset('RelTol',reltol,'AbsTol',abstol);
[T45,X45]=ode45(@fJacVanDerPol,tspan,x0,options,mu);

[T15,X15]=ode15s(@fJacVanDerPol,tspan,x0,options,mu);
disp('Done')

figure()
plot(T, X(:,1), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(T45, X45(:,1), 'DisplayName','ode45')
plot(T15, X15(:,1), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 3)")
xlabel('T')
ylabel('X1')
legend(); hold off

figure()
plot(T, X(:,2), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(T45, X45(:,2), 'DisplayName', 'ode45')
plot(T15, X15(:,2), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 3)")
xlabel('T')
ylabel('X2')
legend(); hold off

figure()
plot(X(:,1), X(:,2), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(X45(:,1), X45(:,2), 'DisplayName', 'ode45')
plot(X15(:,1), X15(:,2), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 3)")
xlabel('X1')
ylabel('X2')
legend(); hold off
%shg % Show graph window
exportgraphics(gcf, append(plotpath, '2_4a.pdf'))

%% 4b) Adaptive on Jac Van Der Pol
tspan = [0, 40];
mu = 20; % mu = 20
x0 = [1.0; 1.0];
h0 = 1/100; % Initial step size
abstol = 1e-06;
reltol = 1e-06;
[T, X] = ExplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

%%% 4b compare with ode45 and ode15
options = odeset('RelTol',reltol,'AbsTol',abstol);
[T45,X45]=ode45(@fJacVanDerPol,tspan,x0,options,mu);

[T15,X15]=ode15s(@fJacVanDerPol,tspan,x0,options,mu);
disp('Done')

figure()
plot(T, X(:,1), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(T45, X45(:,1), 'DisplayName','ode45')
plot(T15, X15(:,1), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 20)")
xlabel('T')
ylabel('X1')
legend(); hold off

figure()
plot(T, X(:,2), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(T45, X45(:,2), 'DisplayName', 'ode45')
plot(T15, X15(:,2), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 20)")
xlabel('T')
ylabel('X2')
legend(); hold off

figure()
plot(X(:,1), X(:,2), 'DisplayName', 'Adaptive Explicit Euler') ; hold on
plot(X45(:,1), X45(:,2), 'DisplayName', 'ode45')
plot(X15(:,1), X15(:,2), 'DisplayName', 'ode15s')
title("Van Der Pol Solution (\mu = 20)")
xlabel('X1')
ylabel('X2')
legend(); hold off
%shg % Show graph window
exportgraphics(gcf, append(plotpath, '2_4b.pdf'))

% CONCLUSION: All three methods are quite similar. No difference can be
% seen between ode15s and ode45, but the (Adaptive) Explicit Euler is
% slightly delayed in the peaks compared to the other two.