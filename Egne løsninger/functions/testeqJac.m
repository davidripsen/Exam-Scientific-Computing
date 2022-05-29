function [dfdx,dgdx]=testeqJac(t,x,lambda)

dfdx = [lambda];
dgdx = eye(1);