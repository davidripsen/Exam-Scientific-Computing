%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 6 - Dormand-Prince 5(4)
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

%% 1. Descripe DOPRI54 with adaptive step size
%%%%% Opg. 3 - Dormand-Prince 5(4) - DOPRI54
% 3.1) Description: Er en Runge-Kutta metode af orden 5,
% dvs. den numeriske approximations fejl (fra den eksakte løsning)
% har orden 5, altså h^0.5, dvs. for lille step size (altid mindre end én),
% bliver fejlen ekstrem lille.

% DOPRI54 bruger en orden 4 'embedded method'. Differencen på løsningen
% mellem orden 5 og orden 4 metoden bruges som estimat på metodens fejl.
% Dette bruges til adaptive stepsize: Op- og nedjustér step sizen h,
% relativt til hvor store fejl, vi tror vi har gang i med vores metode.

%% 2. Implement
solver = ERKSolverErrorEstimationParameters('DOPRI54');

%% 3. On test equation
lambda = -1;
tspan = [0, 32];
x0 = 1;
h = 1/100;
abstol = 1e-05;
reltol = abstol;
solver = ERKSolverErrorEstimationParameters('DOPRI54');
[Tout,Xout,Eout] = AdaptiveERKSolverErrorEstimation(@testeq,tspan,x0,h, ...
    solver,abstol,reltol, lambda);

plot(Tout, Xout)
hold on
plot(Tout, exp(lambda*Tout))
legend()
hold off

%% 4. Van Der Pol
solver = ERKSolverErrorEstimationParameters('DOPRI54');
mus = [3 20];
abstols = [1e-02 1e-06 1e-09];
reltols = abstols;
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
        
        [T,X,E] = AdaptiveERKSolverErrorEstimation(@VanderPolFun,tspan,x0,h0, ...
                solver,abstol,reltol,mu);

        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('Adaptive DOPRI54, tol = %.5g',abstol)) ; hold on
        %shg % Show graph window
    end
        %%% 4b compare with ode45 and ode15
    
    options = odeset('RelTol',reltol,'AbsTol',abstol);
    [T45,X45]=ode45(@VanderPolFun,tspan,x0,options,mu);
    [T15,X15]=ode15s(@VanderPolFun,tspan,x0,options,mu);
    plot(X45(:,1), X45(:,2), DisplayName=sprintf('ode45, tol = %.5g',abstol)')
    plot(X15(:,1), X15(:,2), DisplayName=sprintf('ode15s, tol = %.5g',abstol))
    title(sprintf("Van Der Pol Solution (µ = %i)", mu))
    xlabel('X1')
    ylabel('X2')
    legend('location', 'northwest'); hold off
end
exportgraphics(gcf, append(plotpath, '6_d.pdf'))


%% 5. Adiabatic CSTR








