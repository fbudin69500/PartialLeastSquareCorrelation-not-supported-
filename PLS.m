function [probaInertia,probaSingularValue,Percent,UOutputnames,VOutputnames,U,S,V,Lx,Ly]=PLS(Data,Groups,GroupNorm,splitIndex,nbPerm,names,ProbaPerVector)
    %Normalizing data
    [X,Y]=PLSNormalizeData(Data,GroupNorm,splitIndex);
    %PLS
    R = Y'*X;
    [U,S,V]=svd(R);
    %Permutations
    if ProbaPerVector ~= 0
      [probaInertia,probaSingularValue,Percent]=PLSPermutationsProbaPerVector(X,Y,S,nbPerm);
    else
      [probaInertia,probaSingularValue,Percent]=PLSPermutations(X,Y,S,nbPerm);
    end
    %Find singular values that are below a threshold (p<0.05)
    sI=find(probaSingularValue < 0.05  );
    sI(:,5:end)=[];
    %Boostrap
    [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,GroupNorm,splitIndex,U,V,S,nbPerm,sI);
    UInputnames=names(splitIndex+1:end);
    VInputnames=names(1:splitIndex);
    UOutputnames=cell(size(UInputnames,1),size(sI,2));
    VOutputnames=cell(size(VInputnames,1),size(sI,2));
    for i=1:size(sI,2)
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
    PLSPlot(Lx,size(sI,2),Groups, colors,'Lx=X*V');
    PLSPlot(Ly,size(sI,2),Groups, colors,'Ly=Y*U');
end