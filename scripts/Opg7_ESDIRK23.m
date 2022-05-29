%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 7 - ESDIRK23
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';
fig.Position(3:4)=[550,400]*1.7;

%% 1.
% Calculations + Maple

%% 2. Stability region
% A-stable: Yes
% L-stabke: Yes - ref formler

% myESDIRK23
gamma = 1-1/sqrt(2);
a21 = 1 - sqrt(2)/2;
b1 = sqrt(2)/4;
b2 = b1;
AT = [0 a21 b1;0 gamma b2;0 0 gamma];
c  = [0; 2 - sqrt(2); 1];
b  = AT(:,3);
bhat = [   -sqrt(2)/12 + 1/3; ...
   sqrt(2)/4 + 1/3; ...
    -sqrt(2)/6 + 1/3    ];
d  = b-bhat;
p  = 2;
phat = 3;
s = 3;


alpha = -15:0.02:15;
beta =  -15:0.02:15;
A = AT';

nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);
absR = NaN(nimag, nreal);

for kreal = 1:nreal
   for kimag = 1:nimag
        z = alpha(kreal) + 1i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        f = exp(z);
        absR(kimag,kreal) = abs(R);
   end
end

% ESDIRK23
figure()
fs = 14;
    imagesc(alpha,beta,absR,[0 1]);
    grid on
    colorbar
    axis image
    axis xy
    xlabel('real','fontsize',fs);
    ylabel('imag','fontsize',fs);
    title('|R(z)| - ESDIRK23','fontsize',fs)
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)    
exportgraphics(gcf, append(plotpath, '7_2a.pdf'))

%% 3. Implementation of adaptive ESDIRK23

%% 4. Test on the Van der Pol (and compare with ode45 and ode 15)
Method = "myESDIRK23"

mus = [3 20];
abstols = [1e-02 1e-04 1e-06];
reltols = abstols;
figure('Position', [100, 100, 1200, 600]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);
        %[T, X] = ExplicitEulerAdaptiveStep(...
        %    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        
        [T, X, Gout,info,stats] = ESDIRK(@VanderPolFun,@VanderPolJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,mu);
        
        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('ESDIRK23, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T), info.nFun)) ; hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    sol45=ode45(@VanderPolFun,tspan,x0,options,mu);
    T45 = sol45.x;
    X45 = sol45.y;
    fevals45 = sol45.stats.nfevals;
    sol15=ode15s(@VanderPolFun,tspan,x0,options,mu);
    T15 = sol15.x;
    X15 = sol15.y;
    fevals15 = sol15.stats.nfevals;
    plot(X45(:,1), X45(:,2), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T45), fevals45))
    plot(X15(:,1), X15(:,2), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T15), fevals15))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('X1')
    ylabel('X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '7_4b.pdf'))






%% 7.5 Comparison with other solvers - phase portrait
mus = [3 20];
abstols = [1e-04];
reltols = abstols;
figure('Position', [100, 100, 1200, 600]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);
        subplot(1,2,i)

        % Explicit Euler
        [T, X] = ExplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(X(:,1), X(:,2), DisplayName=sprintf('Explicit Euler Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, length(T)-1)) ; hold on

        % Implicit Euler
        [T,X,iter] = ImplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(X(:,1), X(:,2), DisplayName=sprintf('Implicit Euler Adaptive, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on
        
        % Classical Runge-Kutta
        [T, X, H,fevals] = ClassicalRungeKuttaAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(X(:,1), X(:,2), DisplayName=sprintf('Classical RK Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, fevals)) ; hold on

        % DOPRI54
        solver = ERKSolverErrorEstimationParameters('DOPRI54');
        [T,X,E] = AdaptiveERKSolverErrorEstimation(@VanderPolFun,tspan,x0,h0, ...
                solver,abstol,reltol,mu);
        plot(X(:,1), X(:,2), DisplayName=sprintf('DOPRI54, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on

        % ESDIRK23
        Method = "myESDIRK23"
        [T, X, Gout,info,stats] = ESDIRK(@VanderPolFun,@VanderPolJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,mu);
        plot(X(:,1), X(:,2), DisplayName=sprintf('ESDIRK23, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, info.nFun)) ; hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    sol45=ode45(@VanderPolFun,tspan,x0,options,mu);
    T45 = sol45.x;
    X45 = sol45.y;
    fevals45 = sol45.stats.nfevals;
    sol15=ode15s(@VanderPolFun,tspan,x0,options,mu);
    T15 = sol15.x;
    X15 = sol15.y;
    fevals15 = sol15.stats.nfevals;
    plot(X45(:,1), X45(:,2), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T45), fevals45))
    plot(X15(:,1), X15(:,2), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T15), fevals15))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('X1')
    ylabel('X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '7_5b.pdf'))

%% 7.5 - for X1

mus = [3 20];
abstols = [1e-04];
reltols = abstols;
figure('Position', [100, 100, 1200, 600]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);
        subplot(1,2,i)

        % Explicit Euler
        [T, X] = ExplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,1), DisplayName=sprintf('Explicit Euler Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, length(T)-1)) ; hold on

        % Implicit Euler
        [T,X,iter] = ImplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,1), DisplayName=sprintf('Implicit Euler Adaptive, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on
        
        % Classical Runge-Kutta
        [T, X, H, fevals] = ClassicalRungeKuttaAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,1), DisplayName=sprintf('Classical RK Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, fevals)) ; hold on


        % DOPRI54
        solver = ERKSolverErrorEstimationParameters('DOPRI54');
        [T,X,E] = AdaptiveERKSolverErrorEstimation(@VanderPolFun,tspan,x0,h0, ...
                solver,abstol,reltol,mu);
        plot(T, X(:,1), DisplayName=sprintf('DOPRI54, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on

        % ESDIRK23
        Method = "myESDIRK23"
        [T, X, Gout,info,stats] = ESDIRK(@VanderPolFun,@VanderPolJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,mu);
        plot(T, X(:,1), DisplayName=sprintf('ESDIRK23, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, info.nFun)) ; hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    sol45=ode45(@VanderPolFun,tspan,x0,options,mu);
    T45 = sol45.x;
    X45 = sol45.y;
    fevals45 = sol45.stats.nfevals;
    sol15=ode15s(@VanderPolFun,tspan,x0,options,mu);
    T15 = sol15.x;
    X15 = sol15.y;
    fevals15 = sol15.stats.nfevals;
    plot(T45, X45(1,:), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T45), fevals45))
    plot(T15, X15(1,:), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T15), fevals15))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('T')
    ylabel('X1')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '7_5_X1.pdf'))







%% 7.5 - for X2

mus = [3 20];
abstols = [1e-04];
reltols = abstols;
figure('Position', [100, 100, 1200, 600]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);
        subplot(1,2,i)

        % Explicit Euler
        [T, X] = ExplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,2), DisplayName=sprintf('Explicit Euler Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, length(T)-1)) ; hold on

        % Implicit Euler
        [T,X,iter] = ImplicitEulerAdaptiveStep(...
            @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,2), DisplayName=sprintf('Implicit Euler Adaptive, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on
        
        % Classical Runge-Kutta
        [T, X, H, fevals] = ClassicalRungeKuttaAdaptiveStep(...
            @VanderPolFun,tspan,x0,h0,abstol,reltol,mu);
        plot(T, X(:,2), DisplayName=sprintf('Classical RK Adaptive, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, fevals)) ; hold on

        % DOPRI54
        solver = ERKSolverErrorEstimationParameters('DOPRI54');
        [T,X,E] = AdaptiveERKSolverErrorEstimation(@VanderPolFun,tspan,x0,h0, ...
                solver,abstol,reltol,mu);
        plot(T, X(:,2), DisplayName=sprintf('DOPRI54, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on

        % ESDIRK23
        Method = "myESDIRK23"
        [T, X, Gout,info,stats] = ESDIRK(@VanderPolFun,@VanderPolJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,mu);
        plot(T, X(:,2), DisplayName=sprintf('ESDIRK23, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T)-1, info.nFun)) ; hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    sol45=ode45(@VanderPolFun,tspan,x0,options,mu);
    T45 = sol45.x;
    X45 = sol45.y;
    fevals45 = sol45.stats.nfevals;
    sol15=ode15s(@VanderPolFun,tspan,x0,options,mu);
    T15 = sol15.x;
    X15 = sol15.y;
    fevals15 = sol15.stats.nfevals;
    plot(T45, X45(2,:), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T45), fevals45))
    plot(T15, X15(2,:), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i, fevals = %.i)',abstol, length(T15), fevals15))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('T')
    ylabel('X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '7_5_X2.pdf'))




