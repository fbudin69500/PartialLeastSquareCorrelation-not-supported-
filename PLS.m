function [probaInertia,probaSingularValue,Percent,UOutputnames,VOutputnames,U,S,V,Lx,Ly]=PLS(Data,Groups,GroupNorm,splitIndex,nbPerm,names,ProbaPerVector)
    %Normalizing data
    [X,Y]=PLSNormalizeData(Data,GroupNorm,splitIndex);
    
    
    %Stack groups
    nbgroups=max(Groups);
    for i=1:nbgroups
        index=find( Groups == i);
        data=[X Y];
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
    %PLS
    %R = Y'*X;
    [U,S,V]=svd(R);
    %Permutations
    if ProbaPerVector ~= 0
      [probaInertia,probaSingularValue,Percent]=PLSPermutationsProbaPerVector(X,Y,S,nbPerm,U,V,Groups,splitIndex);
    else
      [probaInertia,probaSingularValue,Percent]=PLSPermutations(X,Y,S,nbPerm,U,V);
    end
    %Find singular values that are below a threshold (p<0.05)
    sI=find(probaSingularValue < 0.05  );
    sIPlot=sI(:,1:5);
    %sI(:,5:end)=[];
    %Boostrap
    UOutputnames=0;
    VOutputnames=0;
%     [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,GroupNorm,splitIndex,U,V,S,nbPerm,sIPlot);
%     UInputnames=names(splitIndex+1:end);
%     VInputnames=names(1:splitIndex);
%     UOutputnames=cell(size(UInputnames,1),size(sIPlot,2));
%     VOutputnames=cell(size(VInputnames,1),size(sIPlot,2));
%     for i=1:size(sIPlot,2)
%         Unames=PLSStable(Uratio,UConfInf,UConfSup,UInputnames,i);
%         for j =1:size(Unames,1)
%           UOutputnames{j,i} = Unames{j,1};
%         end
%         Vnames=PLSStable(Vratio,VConfInf,VConfSup,VInputnames,i);
%         for j =1:size(Vnames,1)
%           VOutputnames{j,i} = Vnames{j,1};
%         end
%     end
    %Latent variables & plots
    Lx=X*V;

    for i=1:nbgroups
        index=find( Groups == i);
        Yg=data(index,splitIndex+1:end);
        Uindex=1:90;
        Uindex=Uindex*i;
        Ug=U(Uindex,:);
        Lyg=Yg*Ug;
        if exist('Ly','var')
            Ly=[Ly;Lyg];
        else
            Ly=Lyg;
        end
    end
    
    
    colors={'bs','bo','mo','ms','gd','c<'};
    %LxPlot=Lx(:,1:5);
    %LyPlot=Ly(:,1:5);
    LxPlot=Lx(:,sIPlot);
    LyPlot=Ly(:,sIPlot);
    PLSPlot(LxPlot,size(LxPlot,2),Groups, colors,'Lx=X*V');
    PLSPlot(LyPlot,size(LyPlot,2),Groups, colors,'Ly=Y*U');
    PLSPlotXY(LxPlot,LyPlot,size(LxPlot,2),Groups, colors,'Lx * Ly');
    PLSPlotXY(Lx(:,1:5),Ly(:,1:5),5,Groups, colors,'Lx * Ly 5 first saliences');
    %Permutations Lx and Ly
    for j=1:nbPerm
        p=randperm(size(LxPlot,1));
        Lxp=Lx(p,1:size(LxPlot,2));
        avgLxp=avg(Lxp);
    end
end