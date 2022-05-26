% Uge8_RungeKutta - Classical Runge Kutta
solver = ERKSolverErrorEstimationParameters('RK44')
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

% a) Describe Classical Runge-Kutta
% A 4-step method wich is a weighed average of in-between steps when
% estimating x_k+1 from x_k.


%% Phase portraint ---- FIXED STEP SIZE
mus = [3 20];
hs = [0.1 0.01 0.001];
figure('Position', [100, 100, 1300, 800]);
for i = 1:length(mus)
    for j = 1:length(hs)
        mu = mus(i);
        h = hs(j);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];

        [T, X] = ClassicalRungeKuttaFixedStep(@VanderPolFun,tspan,x0,h,mu);

        % Choose N = length(T), i.e. same number of steps for the fixed
        % step size Euler.
        %[Tfix, Xfix, fcountfix] = ExplicitEulerFixedStepSize(@VanderPolFun, tspan(1), tspan(2), length(T)-1, x0, mu);
        
        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('Fixed Classical RK, h = %.5g (steps = %.i)',h, length(T)-1)) ; hold on
        %plot(Xfix(:,1), Xfix(:,2), DisplayName=sprintf('Fixed Explicit Euler, h = %.2g (steps = %.i, fevals=%.i)',1/(length(T-1)), length(Tfix), fcountfix)); hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    [T45,X45]=ode45(@VanderPolFun,tspan,x0,options,mu);
    [T15,X15]=ode15s(@VanderPolFun,tspan,x0,options,mu);
    plot(X45(:,1), X45(:,2), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i)',abstol, length(T45)))
    plot(X15(:,1), X15(:,2), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i)',abstol, length(T15)))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlim([min(X(:,1)), max(X(:,1))])
    ylim([min(X(:,2)), max(X(:,2))])
    xlabel('X1')
    ylabel('X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '5_4fix_a.pdf'))






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADAPTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Actually run phase portrait
mus = [3 20];
abstols = [1e-02 1e-04];
reltols = abstols;
figure('Position', [100, 100, 1300, 800]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);

        [T, X, H] = ClassicalRungeKuttaAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);


        % Choose N = length(T), i.e. same number of steps for the fixed
        % step size Euler.
        %[Tfix, Xfix, fcountfix] = ExplicitEulerFixedStepSize(@VanderPolFun, tspan(1), tspan(2), length(T)-1, x0, mu);
        
        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('Adaptive Classical RK, tol = %.5g (steps = %.i)',abstol, length(T))) ; hold on
        %plot(Xfix(:,1), Xfix(:,2), DisplayName=sprintf('Fixed Explicit Euler, h = %.2g (steps = %.i, fevals=%.i)',1/(length(T-1)), length(Tfix), fcountfix)); hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    [T45,X45]=ode45(@VanderPolFun,tspan,x0,options,mu);
    [T15,X15]=ode15s(@VanderPolFun,tspan,x0,options,mu);
    plot(X45(:,1), X45(:,2), DisplayName=sprintf('ode45, tol = %.5g (steps = %.i)',abstol, length(T45)))
    plot(X15(:,1), X15(:,2), DisplayName=sprintf('ode15s, tol = %.5g (steps = %.i)',abstol, length(T15)))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('X1')
    ylabel('X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '5_4.pdf'))

%% Run X1 and X2
mus = [3 20];
abstols = [1e-02 1e-04];
reltols = abstols;
figure('Position', [100, 100, 1300, 800]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);

        [T, X, H] = ClassicalRungeKuttaAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);


        subplot(1,2,i)
        plot(T,X(:,1), DisplayName=sprintf('x_1  (tol = %.5g)',abstol)) ; hold on
        plot(T,X(:,2), DisplayName=sprintf('x_2  (tol = %.5g)',abstol)) ; hold on
        

    end
        %%% 4b compare with ode45 and ode15
    
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('t')
    ylabel('X1 or X2')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '5_4a.pdf'))

%% Run step size H
mus = [3 20];
abstols = [1e-02 1e-04];
reltols = abstols;
figure('Position', [100, 100, 1300, 800]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        %mu = 3; % mu = 20
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);

        [T, X, H] = ClassicalRungeKuttaAdaptiveStep(...
    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);


        subplot(1,2,i)
        plot(T, H, DisplayName=sprintf('Step size  (tol = %.5g)',abstol)) ; hold on

    end    
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('t')
    ylabel('Step size')
    legend('location', 'southoutside'); hold off
end
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '5_4b.pdf'))





