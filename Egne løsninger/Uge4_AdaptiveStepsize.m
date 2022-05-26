%% Uge 4 - Adaptive Step Size

% On Jac Van Der Pol
tspan = [0, 50];
mu = 10;
x0 = [2.0; 0.0];
h0 = 1/100; % Initial step size
abstol = 1e-05;
reltol = 1e-05;
[T, X] = ExplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

h = T(2:end)-T(1:(end-1));
plot(T,X(:,1))
hold on
plot(T,X(:,2))
hold off
figure
%plot(X(:,1), X(:,2), 'ob')
%plot(X(:,1),X(:,2))
%plot(T(2:end)-T(1:(end-1)))
plot(T(2:end),h)


%% Classical RungeKuttaAdaptiveStep (4step)
tspan = [0, 20];
mu = 3;
x0 = [2.0; 0.0];
h0 = 1/10; % Initial step size
abstol = 1e-03;
reltol = 1e-03;
[T, X] = ClassicalRungeKuttaAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

plot(T,X(:,1))
hold on
plot(T,X(:,2))
hold off
%figure
%plot(T(2:end),h)


%% Adaptive Implicit Euler (Exact Newton)

% On Jac Van Der Pol
tspan = [0, 50];
mu = 10;
x0 = [2.0; 0.0];
h0 = 1/100; % Initial step size: 1/100
abstol = 1e-05;
reltol = 1e-05;
[T, X, iter] = ImplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

%h = T(2:end)-T(1:(end-1));
plot(T,X(:,1))
hold on
plot(T,X(:,2))
hold off
figure
plot(T(2:end),h)
fprintf('Number of iterations %i', iter)

%% Inexact Implicit Adaptive Euler

% On Jac Van Der Pol
tspan = [0, 50];
mu = 10;
x0 = [2.0; 0.0];
h0 = 1/100; % Initial step size
abstol = 1e-05;
reltol = 1e-05;
[T, X, iter] = ImplicitEulerAdaptiveInexactNewt(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

h = T(2:end)-T(1:(end-1));
plot(T,X(:,1))
hold on
plot(T,X(:,2))
hold off
figure
plot(T(2:end),h)
plot(X(:,1), X(:,2))
fprintf('Number of iterations %i', iter)

%% Compare solutions (of Exact and Inexact Newton)
[Ta, Xa, itera] = ImplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
[Tb, Xb, iterb] = ImplicitEulerAdaptiveInexactNewt(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

Xa(end,1)
Xb(end,1)

Xa(end,2)
Xb(end,2)
