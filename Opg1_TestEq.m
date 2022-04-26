%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 1 - Test equation for ODEs
clear clc
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

% x'(t) = lambda*x(t), where lambda = -1, x(0) = x0 = 1

%% 1) Analytical solution
% Try x(t) = exp(lambda*t); Then x'(t) = lambda*exp(lambda*t) = lambda*x(t)

%% 2) Local and global truncation error
% Se slides for formler.
% Local: One step error fra analytisk løsning, hvor analytisk løsning
% iterativt flyttes (fejlagtigt) til starten  af næste step

% Global: Samlede propagerede fejl fra analytisk, i sidste step.

%% 3) Computation of local and global truncation errors
% løst i hhv. 4) og 5).

%% 4) Local Error and ORDER
% Order = antal Taylor-led der inkluderes. Local error er proportional med
% step size ^ (order + 1)
x0 = 1;
lambda = -1;
T0 = 0;
Tf = 1;
Ns = round(linspace(10^1, 10^2.6, 100));
hs = (Tf-T0) ./ Ns;
E = zeros(length(hs),3);
tspan = [T0, Tf];

% Calculate local error in 1st step (i.e. X(2))
for i=1:length(Ns)
    [T_exp, X_exp] = ExplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    [T_imp, X_imp] = ImplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    solver = ERKSolverErrorEstimationParameters('RK44');
    [T_cRK, X_cRK,~] = ERKSolverErrorEstimation(@testeq, tspan, x0, hs(i), solver, lambda);
    
    E(i,1) = abs(X_exp(2) - x0*exp(hs(i)*lambda));
    E(i,2) = abs(X_imp(2) - x0*exp(hs(i)*lambda));
    E(i,3) = abs(X_cRK(2) - x0*exp(hs(i)*lambda));
end

loglog(hs, abs(E(:,1)), 'Marker','*', 'DisplayName', 'Explicit Euler');
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

hold off;
legend()
shg % Show current figure
exportgraphics(gcf, append(plotpath, '1_4.pdf'))

%% 5) Global Error
% ========================================================================
x0 = 1;
lambda = -1;
T0 = 0;
Tf = 1;
Ns = round(linspace(10^1, 10^3, 100));
hs = (Tf-T0) ./ Ns;
E = zeros(length(hs),3);
tspan = [T0, Tf];


for i=1:length(Ns)
    [T_exp, X_exp] = ExplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    [T_imp, X_imp] = ImplicitEulerFixedStepSize(@testeq, T0, Tf, Ns(i), x0, lambda);
    solver = ERKSolverErrorEstimationParameters('RK44');
    [T_cRK, X_cRK,~] = ERKSolverErrorEstimation(@testeq, tspan, x0, hs(i), solver, lambda);
    
    E(i,1) = abs(X_exp(end) - x0*exp(Tf*lambda));
    E(i,2) = abs(X_imp(end) - x0*exp(Tf*lambda));
    E(i,3) = abs(X_cRK(end) - x0*exp(Tf*lambda));
end

loglog(hs, abs(E(:,1)), 'Marker','*', 'DisplayName', 'Explicit Euler');
grid on
xlabel('h')
ylabel('Global error')
title('Global error vs step size')
polyfit(log(hs), log(abs(E(:,1))), 1)

hold on;
loglog(hs, abs(E(:,2)), 'Marker','*', 'DisplayName','Implicit Euler');
polyfit(log(hs), log(abs(E(:,2))), 1)

loglog(hs, abs(E(:,3)), 'Marker','*', 'DisplayName','Classical RK');
polyfit(log(hs), log(abs(E(:,3))), 1)

hold off;
legend()
shg % Show graph window
exportgraphics(gcf, append(plotpath, '1_5.pdf'))

%plot(hspan, E)
% I.e. a clear linear line, with slope ≈ 4. I.e. order = 3 (??? mangler) <=> global error
% is proportional to h^4.



%% 1.6 Explain stability of a method
% bla bla
clear; close all;
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

method = 'RK44';
alpha = -5:0.01:5;
beta =  -5:0.01:5;
solver = ERKSolverErrorEstimationParameters(method);
A = solver.AT';
b = solver.b;
c = solver.d;
d = solver.d;
nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);
absR = NaN(nimag, nreal);
absR_exp = NaN(nimag, nreal);
absR_imp = NaN(nimag, nreal);

for kreal = 1:nreal
   for kimag = 1:nimag
        z = alpha(kreal) + i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        R_exp = 1 + z;
        R_imp = 1 / (1-z);
        f = exp(z);
        absR(kimag,kreal) = abs(R);
        absR_exp(kimag,kreal) = abs(R_exp);
        absR_imp(kimag,kreal) = abs(R_imp);
   end
end

% Classical Runge-Kutta
figure()
fs = 14;
    imagesc(alpha,beta,absR,[0 1]);
    grid on
    colorbar
    axis image
    axis xy
    xlabel('real','fontsize',fs);
    ylabel('imag','fontsize',fs);
    title('|R(z)| - Classical Runge-Kutta','fontsize',fs)
exportgraphics(gcf, append(plotpath, '1_6a.pdf'))

% Explicit Euler
figure()
fs = 14;
    imagesc(alpha,beta,absR_exp,[0 1]);
    grid on
    colorbar
    axis image
    axis xy
    xlabel('real','fontsize',fs);
    ylabel('imag','fontsize',fs);
    title('|R(z)| - Explicit Euler','fontsize',fs)
exportgraphics(gcf, append(plotpath, '1_6b.pdf'))

% Implicit Euler
figure()
fs = 14;
    imagesc(alpha,beta,absR_imp,[0 1]);
    grid on
    colorbar
    axis image
    axis xy
    xlabel('real','fontsize',fs);
    ylabel('imag','fontsize',fs);
    title('|R(z)| - Implicit Euler','fontsize',fs)
exportgraphics(gcf, append(plotpath, '1_6c.pdf'))

% Comments:
% A-stable: For all z where R(z) < 1, then Re(z) < 0
% the case for: 

% L-stable: A-stable + For z -> -inf, R(z) = 0
% the case for:











