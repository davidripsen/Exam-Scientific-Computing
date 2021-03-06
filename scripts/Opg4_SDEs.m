%%% 02686 Scientific Computing for Differential Equations - Exam
%%% by David Ribberholt Ipsen (s164522)
%%% Exercise 4 - Solvers for SDEs
plotpath = '/Users/davidipsen/Documents/DTU/4. Semester (MSc)/Scientific Computing/Exam-Scientific-Computing/plots/';

%% 1. Multivariate Standard Wiener
% se StdWienerProcess

%% 2. (explicit-explicit fixed step size) - state dependent diffusion
% Change VanDerPolDiffusion2 -> 1 for state-independent diffusion


mus = [3 20];
sigmas = [0.1 2];

for mu = mus
    j = 0;
    figure()
    for sigma = sigmas
        %%%%%%%%%%%%%%% EXP-EXP
        j = j+1;
        x0 = [1.0; 1.0];
        p = [mu; sigma];
        tf = 32;
        nw = 1;
        N = 32/0.001; % Length of realization
        Ns = 5; % Number of realizations
        seed = 100;
         
        [W,T,~]=StdWienerProcess(tf,N,nw,Ns,seed);
        X = zeros(length(x0),N+1,Ns);
        for i=1:Ns
            X(:,:,i) = SDEsolverExplicitExplicit(...
                @VanderPolDrift, @VanderPolDiffusion1,...
                T,x0,W(:,:,i),p);
        end
        % Deterministic by setting sigma=0
        Xd = SDEsolverExplicitExplicit(...
                @VanderPolDrift,@VanderPolDiffusion1,...
                T,x0,W(:,:,i),[mu; 0.0]);
        % Plot
        subplot(2,2,j)
        hold on
        for realization = 1:Ns
            plot(X(1,:,realization), X(2,:,realization))
        end
        plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
        title(sprintf('Exp-Exp  (µ = %.i, σ = %.1f)', mu, sigma))
        hold off
        %figure

        %%%%%%%%%%%%%% REPEAT FOR IMP-EXP
        j = j+1;
        x0 = [1.0; 1.0];
        p = [mu; sigma];
        tf = 32;
        nw = 1;
        N = 5000; % Length of realization
        Ns = 5; % Number of realizations
        seed = 100;
         
%         [W,T,~]=StdWienerProcess(tf,N,nw,Ns,seed);
%         X = zeros(length(x0),N+1,Ns);
%         for i=1:Ns
%             X(:,:,i) = SDEsolverImplicitExplicit(...
%                 @VanderPolDrift, @VanderPolDiffusion1,...
%                 T,x0,W(:,:,i),p);
%         end
        
        % Deterministic by setting sigma=0
        Xd = SDEsolverImplicitExplicit(...
                @VanderPolDrift,@VanderPolDiffusion1,...
                T,x0,W(:,:,i),[mu; 0.0]);

        % Plot
        subplot(2,2,j)
        hold on
        for realization = 1:Ns
            plot(X(1,:,realization), X(2,:,realization))
        end
        plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
        title(sprintf('Imp-Exp  (µ = %.i, σ = %.1f)', mu, sigma))
        hold off
        %figure
        %set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
    end
end
%exportgraphics(gcf, append(plotpath, '4a.pdf'))

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



%% 3. (implicit-explicit drift-diffusion solve, fixed step size
mus = [3 20];
sigmas = [0.1 2];

for mu = mus
    j = 0;
    figure()
    for sigma = sigmas
        %%%%%%%%%%%%%%% EXP-EXP
        j = j+1;
        x0 = [1.0; 1.0];
        p = [mu; sigma];
        tf = 32;
        nw = 1;
        N = 32/0.001; % Length of realization
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
        subplot(2,2,j)
        hold on
        for realization = 1:Ns
            plot(X(1,:,realization), X(2,:,realization))
        end
        plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
        title(sprintf('Exp-Exp  (µ = %.i, σ = %.1f)', mu, sigma))
        hold off
        %figure

        %%%%%%%%%%%%%% REPEAT FOR IMP-EXP
        j = j+1;
        x0 = [1.0; 1.0];
        p = [mu; sigma];
        tf = 32;
        nw = 1;
        N = 5000; % Length of realization
        Ns = 5; % Number of realizations
        seed = 100;
         
%         [W,T,~]=StdWienerProcess(tf,N,nw,Ns,seed);
%         X = zeros(length(x0),N+1,Ns);
%         for i=1:Ns
%             X(:,:,i) = SDEsolverImplicitExplicit(...
%                 @VanderPolDrift, @VanderPolDiffusion1,...
%                 T,x0,W(:,:,i),p);
%         end
        
        % Deterministic by setting sigma=0
        Xd = SDEsolverImplicitExplicit(...
                @VanderPolDrift,@VanderPolDiffusion2,...
                T,x0,W(:,:,i),[mu; 0.0]);

        % Plot
        subplot(2,2,j)
        hold on
        for realization = 1:Ns
            plot(X(1,:,realization), X(2,:,realization))
        end
        plot(Xd(1,:,1), Xd(2,:,1), 'linewidth', 2)
        title(sprintf('Imp-Exp  (µ = %.i, σ = %.1f)', mu, sigma))
        hold off
        %figure
        %set(findall(0, '-property', 'fontsize'), 'fontsize', 17)
    end
end