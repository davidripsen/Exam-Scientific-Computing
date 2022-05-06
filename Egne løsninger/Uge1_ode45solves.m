%%Driver
mu = 10;
x0 = [2.0; 0.0];
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode45(@VanDerPol,[0 5*mu],x0,options,mu);

%% Plot
figure;
plot(X(:,1), X(:,2))
figure;
plot(T, X(:,1))
figure;
plot(T, X(:,2))


%% Driver - ode15s
mu = 10;
x0 = [2.0; 0.0];
options = odeset('Jacobian',@JacVanDerPol,'RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode15s(@VanDerPol,[0 5*mu],x0,options,mu);
disp('Done')

% And plot
plot(X(:,1), X(:,2))
%plot(T, X(:,1))
%plot(T, X(:,2))

%% Predator-Prey Model
a = 1;
b = 1;
x0 = [2; 2];
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode45(@PreyPredator,[0 50],x0,options,a,b);

% and plot
plot(X(:,1), X(:,2))

%% Lorentz Attractor Model
sigma= 10;
rho= 28;
beta= 8/3;
eta = sqrt(rho*(beta-1));
x0 = [rho-1; eta; eta] + [0; 0; 3];

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode45(@Lorentz,[0 50],x0,options, sigma,rho,beta);

% and plot
figure;
plot(T, X(:,1))
figure;
plot(T, X(:,2))
figure;
plot(T, X(:,3))


%% Functions (in order of appereance)
function xdot = VanDerPol(t,x,mu)
     % VANDERPOL  Implementation of the Van der Pol model
     %
     % Syntax: xdot = VanDerPol(t,x,mu)
     xdot=zeros(2,1);
     xdot(1) = x(2);
     xdot(2) = mu*(1-x(1)*x(1))*x(2)-x(1);
end

function Jac = JacVanDerPol(t,x,mu)
    % JACVANDERPOL  Jacobian for the Van der Pol Equation
    %
    % Syntax: Jac = JacVanDerPol(t,x,mu)
    Jac = zeros(2,2);
    Jac(2,1) = -2*mu*x(1)*x(2)-1.0;
    Jac(1,2) = 1.0;
    Jac(2,2) = mu*(1-x(1)*x(1));
end


function xdot = PreyPredator(t,x,a,b)
    % PREYPREDATOR  The Prey-Predator Model
    %
    % Syntax: xdot = PreyPredator(t,x,a,b)
    xdot = zeros(2,1);
    xdot(1) = a*(1-x(2))*x(1);
    xdot(2) = -b*(1-x(1))*x(2);
end

function xdot = Lorentz(t,x,sigma,rho,beta)
    % Lorentz Attractor model
    xdot = zeros(3,1);
    xdot(1) = sigma * (x(2)-x(1));
    xdot(2) = x(1)*(rho - x(3)) - x(2);
    xdot(3) = x(1)*x(2) - beta*x(3);
end