% Uge 6 - Stochastic Differential Equations

%% First: ODE Adaptive Step Size

% On Jac Van Der Pol
tspan = [0, 50];
mu = 10;
x0 = [2.0; 0.0];
h0 = 1/100; % Initial step size
abstol = 1e-05;
reltol = 1e-05;
[T, X] = ExplicitEulerAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);

%h = T(2:end)-T(1:(end-1));
%plot(T,X(:,1))
%hold on
%plot(T,X(:,2))
%hold off
figure
plot(X(:,1), X(:,2), 'ob')
%plot(X(:,1),X(:,2))
%plot(T(2:end)-T(1:(end-1)))
%plot(T(2:end),h)


%% Now: SDE (explicit-explicit fixed step size) - state dependent diffusion
% Change VanDerPolDiffusion2 -> 1 for state-independent diffusion


mu = 3;
sigma = 1;
x0 = [0.5; 0.5];
p = [mu; sigma];
tf = 5*mu;
nw = 1;
N = 1000; % Length of realization
Ns = 5; % Number of realizations
seed = 100;

[W,T,~]=StdWienerProcess(tf,N,nw,Ns,seed);
X = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X(:,:,i) = SDEsolverExplicitExplicit(...
        @VanderPolDrift, @VanderPolDiffusion2,...
        T,x0,W(:,:,i),p);
end

% Deterministic by setting sigma=0
Xd = SDEsolverExplicitExplicit(...
        @VanderPolDrift,@VanderPolDiffusion2,...
        T,x0,W(:,:,i),[mu; 0.0]);

% Plot
hold on
for realization = 1:Ns
    plot(X(1,:,realization), X(2,:,realization))
end
plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
hold off
%figure

% hold on
% subplot(2,1,1);
% for realization = 1:Ns
%     plot(X(1,:,realization));
% end
% hold off
% 
% hold on
% subplot(2,1,2);
% for realization = 1:Ns
%     plot(X(2,:,realization));
% end
% hold off



%% SDE (implicit-explicit drift-diffusion solve, fixed step size

mu = 3;
sigma = 0.5;
x0 = [0.5; 0.5];
p = [mu; sigma];
tf = 5*mu;
nw = 1;
N = 1000; % Length of realization
Ns = 5; % Number of realizations
seed = 100;

[W,T,~]=StdWienerProcess(tf,N,nw,Ns,seed);
X = zeros(length(x0),N+1,Ns);
for i=1:Ns
    X(:,:,i) = SDEsolverImplicitExplicit(...
        @VanderPolDrift, @VanderPolDiffusion2,...
        T,x0,W(:,:,i),p);
end

% Deterministic by setting sigma=0
Xd = SDEsolverImplicitExplicit(...
        @VanderPolDrift,@VanderPolDiffusion2,...
        T,x0,W(:,:,i),[mu; 0.0]);

% Plot
hold on
for realization = 1:Ns
    plot(X(1,:,realization), X(2,:,realization))
end
plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
hold off
