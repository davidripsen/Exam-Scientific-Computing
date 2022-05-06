% Uge9_RungeKutta+adaptive

%%%%% Opg. 3 - Dormand-Prince 5(4) - DOPRI54
% 3.1) Description: Er en Runge-Kutta metode af orden 5,
% dvs. den numeriske approximations fejl (fra den eksakte løsning)
% har orden 5, altså h^0.5, dvs. for lille step size (altid mindre end én),
% bliver fejlen ekstrem lille.

% DOPRI54 bruger en orden 4 'embedded method'. Differencen på løsningen
% mellem orden 5 og orden 4 metoden bruges som estimat på metodens fejl.
% Dette bruges til adaptive stepsize: Op- og nedjustér step sizen h,
% relativt til hvor store fejl, vi tror vi har gang i med vores metode.

% 3.2) Implement
solver = ERKSolverErrorEstimationParameters('DOPRI54');

% Test equation
lambda = -1;
tspan = [0, 30];
x0 = 1;
h = 1/100;
[Tout,Xout,Eout] = AdaptiveERKSolverErrorEstimation(@testeq,tspan,x0,h, ...
    solver,lambda);

plot(Tout, Xout)
hold on
plot(Tout, exp(lambda*Tout))
legend()
hold off