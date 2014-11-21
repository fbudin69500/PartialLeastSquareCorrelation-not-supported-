function [probaInertia,probaSingularValue,Percent,UOutputnames,VOutputnames]=PLS(Data,Groups,GroupNorm,splitIndex,nbPerm,names)
    %Normalizing data
    [X,Y]=PLSNormalizeData(Data,GroupNorm,splitIndex);
    %PLS
    R = Y'*X;
    [U,S,V]=svd(R);
    %Permutations
    [probaInertia,probaSingularValue,Percent]=PLSPermutations(X,Y,S,nbPerm);
    %Find singular values that are below a threshold (p<0.05)
    sI=find(probaSingularValue < 0.05  );
    maxSI=max(sI);
    %Boostrap
    [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,GroupNorm,splitIndex,U,V,S,nbPerm,maxSI);
    UInputnames=names(splitIndex+1:end);
    VInputnames=names(1:splitIndex);
    UOutputnames=cell(size(UInputnames,1),maxSI);
    VOutputnames=cell(size(VInputnames,1),maxSI);
    for i=1:maxSI
        Unames=PLSStable(Uratio,UConfInf,UConfSup,UInputnames,i);
        for j =1:size(Unames,1)
          UOutputnames{j,i} = Unames{j,1};
        end
        Vnames=PLSStable(Vratio,VConfInf,VConfSup,VInputnames,i);
        for j =1:size(Vnames,1)
          VOutputnames{j,i} = Vnames{j,1};
        end
    end
    %Latent variables & plots
    Lx=X*V;
    Ly=Y*U;
    colors={'bs','bo','mo','ms','gd','c<'};
    PLSPlot(Lx,maxSI,Groups, colors,'Lx=X*V');
    PLSPlot(Ly,maxSI,Groups, colors,'Ly=Y*U');
end