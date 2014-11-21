% Compute permutation: sampling without replacement
function [I,S,Xp]=MyPermute(X,Y)
sizeX=size(X);
p=randperm(sizeX(1));
Xp=X(p,1:sizeX(2));
[U,S,V]=svd(Y'*Xp);
I=sum(sum(S));
end