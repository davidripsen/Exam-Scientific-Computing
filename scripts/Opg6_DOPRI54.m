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
tspan = [0, 10];
x0 = 1;
h = 1/100;

figure('Position', [100, 100, 1300, 800]);
subplot(2,1,1)
tols = [1e-02 1e-04 1e-06];
for i=1:length(tols)
    abstol = tols(i);
    reltol = abstol;
    solver = ERKSolverErrorEstimationParameters('DOPRI54');
    [Tout,Xout,Eout] = AdaptiveERKSolverErrorEstimation(@testeq,tspan,x0,h, ...
        solver,abstol,reltol, lambda);
    
    plot(Tout, Xout, 'DisplayName', sprintf('DOPRI54  (tol=%.2g)', abstol))
    hold on
    disp(Tout)
end
plot(linspace(0,10,1000), exp(lambda*linspace(0,10,1000)), 'DisplayName','Exact')
title("Test equation (lambda=1)")
xlabel("t")
ylabel("x")
legend()
hold off
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)

%%% 3. On test equation - 2
subplot(2,1,2)

tols = [1e-02 1e-04 1e-06];
for i=1:length(tols)
    abstol = tols(i);
    reltol = abstol;
    solver = ERKSolverErrorEstimationParameters('DOPRI54');
    [Tout,Xout,Eout] = AdaptiveERKSolverErrorEstimation(@testeq,tspan,x0,h, ...
        solver,abstol,reltol, lambda);
  
    plot(Tout, abs(exp(lambda*Tout) - Xout) - Eout, 'DisplayName', sprintf('DOPRI54  (tol=%.2g)', abstol))
    hold on
end
%xlim([0 0.1])
%ylim([0 1])
title("Test equation (lambda=1)")
xlabel("t")
ylabel("e - e_{estimated}")
legend()
hold off
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
exportgraphics(gcf, append(plotpath, '6_3a.pdf'))

%% 3b. Order
% Local Error and ORDER
% Order = antal Taylor-led der inkluderes. Local error er proportional med
% step size ^ (order + 1)
x0 = 1;
lambda = -1;
T0 = 0;
Tf = 1;
Ns = round(linspace(10^1, 10^2, 100));
hs = (Tf-T0) ./ Ns;
E = zeros(length(hs),4);
tspan = [T0, Tf];

% Calculate local error in 1st step (i.e. X(2))
for i=1:length(Ns)
    [T_exp, X_exp] = ExplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    [T_imp, X_imp] = ImplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    solver = ERKSolverErrorEstimationParameters('RK44');
    [T_cRK, X_cRK,~] = ERKSolverErrorEstimation(@testeq, tspan, x0, hs(i), solver, lambda);

    solver = ERKSolverErrorEstimationParameters('DOPRI54');
    [Tout,Xout,Eout] = ERKSolverErrorEstimation(@testeq, tspan, x0, hs(i), solver, lambda)
    
    E(i,1) = abs(X_exp(2) - x0*exp(hs(i)*lambda));
    E(i,2) = abs(X_imp(2) - x0*exp(hs(i)*lambda));
    E(i,3) = abs(X_cRK(2) - x0*exp(hs(i)*lambda));
    
    E(i,4) = abs(Xout(2) - x0*exp(hs(i)*lambda));
    disp( (Tout(2)-Tout(1)) - hs(i))
end

loglog(hs, abs(E(:,1)), 'Marker','*', 'DisplayName', 'Explicit Euler');
ylim([1e-16 1e-00])
grid on
xlabel('h')
ylabel('Local error')
title('Local error vs step size')
polyfit(log(hs), log(abs(E(:,1))), 1)

hold on;
loglog(hs, abs(E(:,2)), 'Marker','*', 'DisplayName','Implicit Euler');
polyfit(log(hs), log(abs(E(:,2))), 1)

loglog(hs, abs(E(:,3)), 'Marker','*', 'DisplayName','Classical RK');
polyfit(log(hs), log(abs(E(:,3))), 1)

loglog(hs, abs(E(:,4)), 'Marker','*', 'DisplayName','DOPRI54');
polyfit(log(hs), log(abs(E(:,4))), 1)

hold off;
legend('Location','northwest')
shg % Show current figure
exportgraphics(gcf, append(plotpath, '6_3b.pdf'))
% Other methods for reference.
% slope (DOPRI54) = 6.03 <=> local error is proportional to (step size)^6
% (which is a really low number for step size < 1) <=> order = 5

%% 3c) Stability
% bla bla
clear; close all;
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

alpha = -5:0.01:5;
beta =  -5:0.01:5;

solver_dopri = ERKSolverErrorEstimationParameters('DOPRI54');
d_dopri = solver_dopri.d;
A_dopri = solver_dopri.AT';
b_dopri = solver_dopri.b;

nreal = length(alpha);
nimag = length(beta);

I = eye(size(A));
e = ones(size(A,1),1);

I_dopri = eye(size(A_dopri));
e_dopri = ones(size(A_dopri,1),1);

absR = NaN(nimag, nreal);
absR_exp = NaN(nimag, nreal);
absR_imp = NaN(nimag, nreal);
absR_dopri = NaN(nimag, nreal);

absR = NaN(nimag, nreal);
absEhatmE = NaN(nimag, nreal);
absEhat   = NaN(nimag, nreal);
absE = NaN(nimag, nreal);
absF = NaN(nimag, nreal);

for kreal = 1:nreal
   for kimag = 1:nimag
        z = alpha(kreal) + i*beta(kimag);
        tmp_dopri = (I_dopri-z*A_dopri)\e_dopri;
        R = 1 + z*b_dopri'*tmp_dopri;
        Ehat = z*d_dopri'*tmp_dopri;
        f = exp(z);
        E = R-f;
        EhatmE = Ehat-E;
        absR(kimag,kreal) = abs(R);
        absEhatmE(kimag,kreal) = abs(EhatmE);
        absEhat(kimag,kreal)   = abs(Ehat);
        absE(kimag,kreal) = abs(E);
        absF(kimag,kreal) = abs(f);
   end
end

fs = 14;
    imagesc(alpha,beta,absR,[0 1]);
    %imagesc(alpha,beta,absEhatmE,[0 1]);
    %imagesc(alpha,beta,absEhat,[0 1]);
    %imagesc(alpha,beta,absE,[0 1]);
    grid on
    colorbar
    axis image
    axis xy
    xlabel('real','fontsize',fs);
    ylabel('imag','fontsize',fs);
    title('|R(z)| - DOPRI54','fontsize',fs)
exportgraphics(gcf, append(plotpath, '6_3c.pdf'))

% Comments:
% A-stable: For all z where R(z) < 1, then Re(z) < 0
% the case for: 

% L-stable: A-stable + For z -> -inf, R(z) = 0
% the case for:




%% 4. Van Der Pol
solver = ERKSolverErrorEstimationParameters('DOPRI54');
mus = [3 20];
abstols = [1e-02 1e-04 1e-06];
reltols = abstols;
figure('Position', [100, 100, 1300, 800]);
for i = 1:length(mus)
    for j = 1:length(abstols)
        mu = mus(i);
        tspan = [0, 32];
        x0 = [1.0; 1.0];
        h0 = 1/100; % Initial step size
        abstol = abstols(j);
        reltol = reltols(j);
        %[T, X] = ExplicitEulerAdaptiveStep(...
        %    @fJacVanDerPol,tspan,x0,h0,abstol,reltol,mu);
        
        [T,X,E] = AdaptiveERKSolverErrorEstimation(@VanderPolFun,tspan,x0,h0, ...
                solver,abstol,reltol,mu);

        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('Adaptive DOPRI54, tol = %.5g (steps = %.i)',abstol, length(T)-1)) ; hold on
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
    legend('location', 'southoutside'); hold off
end
exportgraphics(gcf, append(plotpath, '6_5.pdf'))
set(findall(0, '-property', 'fontsize'), 'fontsize', 17)

%% 5. Adiabatic CSTR








