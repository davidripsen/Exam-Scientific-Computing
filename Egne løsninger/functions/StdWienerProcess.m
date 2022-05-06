function [W,Tw,dW] = StdWienerProcess(T,N,nW,Ns,seed)
%Standard Wiener process in [0,T]
%Time points
%White noise used to generate the Wiener process
%Final time
%Number of intervals
%Dimension of W(k)
%Number of realizations
%To set the random number generator (optional)

if nargin == 4
    rng(seed);
end
dt = T/N;
dW = sqrt(dt)*randn(nW,N,Ns);
W = [zeros(nW,1,Ns) cumsum(dW,2)];
Tw = 0:dt:T;