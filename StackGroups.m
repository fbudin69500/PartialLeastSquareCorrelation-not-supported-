function [R]=StackGroups(X,Y,Groups)
nbgroups=max(Groups);
for i=1:nbgroups
    index=find( Groups == i);
    Xg=X(index,:);
    Yg=Y(index,:);
    Rg = Yg'*Xg;
    if exist('R','var')
        R=[R;Rg];
    else
        R=Rg;
    end
end