%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 7 - ESDIRK23
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

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
        z = alpha(kreal) + i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        f = exp(z);
        absR(kimag,kreal) = abs(R);
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
    title('|R(z)| - ESDIRK23','fontsize',fs)
exportgraphics(gcf, append(plotpath, '7_2a.pdf'))

%% 3. Implementation of adaptive ESDIRK23

%% 4. Test on the Van der Pol (and compare with ode45 and ode 15)
Method = "myESDIRK23"

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
        
        [T, X, Gout,info,stats] = ESDIRK(@VanderPolFun,@VanderPolJac,tspan(1),tspan(end),x0,h0,abstol,reltol,Method,mu);
        
        subplot(1,2,i)
        plot(X(:,1), X(:,2), DisplayName=sprintf('Adaptive ESDIRK23, tol = %.5g',abstol)) ; hold on
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
exportgraphics(gcf, append(plotpath, '7_4b.pdf'))
