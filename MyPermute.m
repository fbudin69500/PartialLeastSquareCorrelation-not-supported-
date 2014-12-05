% Compute permutation: sampling without replacement
function [I,S,Xp]=MyPermute(X,Y,Uo,Vo,So)
sizeX=size(X);
p=randperm(sizeX(1));
Xp=X(p,1:sizeX(2));
M=Y'*Xp;
[U2,S2,V2]=svd(M'*[Uo*So*Vo]);
R=U2*V2';
[U,S,V]=svd(M*R);
I=sum(sum(S));
end