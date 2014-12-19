% Compute permutation: sampling without replacement
function [I,S,Xp]=MyPermute(X,Y,Groups,splitIndex)
sizeX=size(X);
p=randperm(sizeX(1));
Xp=X(p,1:sizeX(2));

nbgroups=max(Groups);
for i=1:nbgroups
        index=find( Groups == i);
        data=[Xp Y];
        M=data(index,:);
        Xg=M(:,1:splitIndex);
        Yg=M(:,splitIndex+1:end);
        Rg = Yg'*Xg;
        if exist('R','var')
            R=[R;Rg];
        else
            R=Rg;
        end
    end

[U,S,V]=svd(R);
I=sum(sum(S));
end