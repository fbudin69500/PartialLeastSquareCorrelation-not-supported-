% Compute permutation: sampling without replacement
function [I,S,Xp]=MyPermute(X,Y,Groups)
sizeX=size(X);
p=randperm(sizeX(1));
Xp=X(p,1:sizeX(2));

R=StackGroups(Xp,Y,Groups);
try
  [U,S,V]=svd(R);
  I=sum(sum(S));
catch err
  I=-2;
  S=0;
  Xp=0;
  disp('warning, skipped one permutation')
end

end