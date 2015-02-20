function [probaInertia,probaSingularValue,Percent,UOutputnames,VOutputnames,U,S,V,Lx,Ly]=PLS(Data,Groups,GroupNorm,splitIndex,nbPerm,names,ProbaPerVector,useSI)
    %Normalizing data
    nbgroups=max(Groups);
    normedData=PLSNormalizeData(Data,GroupNorm);
    %Split data
    X=normedData(:,1:splitIndex);
    Y=normedData(:,splitIndex+1:end);
    UInputnames=names(splitIndex+1:end);
    VInputnames=names(1:splitIndex);
    %Stack groups
    R=StackGroups(X,Y,Groups);
    %PLS
    %R = Y'*X;
    [U,S,V]=svd(R,'econ');
    %Permutations
    if ProbaPerVector ~= 0
      [probaInertia,probaSingularValue,Percent]=PLSPermutationsProbaPerVector(X,Y,S,nbPerm,Groups);
    else
      [probaInertia,probaSingularValue,Percent]=PLSPermutations(X,Y,S,nbPerm,Groups);
    end
    %Find singular values that are below a threshold (p<0.05)
    if useSI ~= 0
      sI=find(probaSingularValue < 0.05  );
    else
      sI=1:5
    end
    nbSignificant=min(5,size(sI,2));
    sIPlot=sI(:,1:nbSignificant);
    %Latent variables & plots
    Lx=X*V;

    for i=1:nbgroups
        index=find( Groups == i);
        Yg=normedData(index,splitIndex+1:end);
        Uindex=1:size(Y,2);
        Uindex=Uindex+size(Y,2)*(i-1);
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
    PLSPlotXY(Lx(:,1:nbSignificant),Ly(:,1:nbSignificant),nbSignificant,Groups, colors,'Lx * Ly first saliences');
    %Plot LV
    figure('Name','LV - V');
    hold on
    for n=1:size(sIPlot,2)
        subplot(size(sIPlot,2),1,n );
        bar(V(:,sIPlot(n)));
        title([' V ' num2str(sIPlot(n))]);
    end
    hold off
    figure('Name','LV - U');
    hold on
    for n=1:size(sIPlot,2)
        for i=1:nbgroups
            Uindex=1:size(Y,2);
            Uindex=Uindex+size(Y,2)*(i-1);
            Uext=U(Uindex,sIPlot(n));
            if exist('Ureshaped','var')
                Ureshaped=[Ureshaped,Uext];
            else
                Ureshaped=Uext;
            end
        end
        subplot(size(sIPlot,2),1,n );
        bar(Ureshaped);
        title([' U ' num2str(sIPlot(n))]);
        clear Ureshaped;
    end
    hold off

    
    %Boostrap
    UOutputnames=0;
    VOutputnames=0;
    [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,GroupNorm,splitIndex,U,V,S,nbPerm,sIPlot);
    
    UOutputnames=cell(size(UInputnames,1),size(sIPlot,2)*nbgroups);
    VOutputnames=cell(size(VInputnames,1),size(sIPlot,2));
    for i=1:size(sIPlot,2)
        for k=1:nbgroups
            Uindex=1:size(Y,2);
            Uindex=Uindex+(k-1)*size(Y,2);
            Unames=PLSStable(Uratio(Uindex,:),UConfInf(Uindex,:),UConfSup(Uindex,:),UInputnames,i);
            for j =1:size(Unames,1)
               UOutputnames{j,i*(nbgroups-1)+(k-1)} = Unames{j,1};
            end
        end
        Vnames=PLSStable(Vratio,VConfInf,VConfSup,VInputnames,i);
        for j =1:size(Vnames,1)
          VOutputnames{j,i} = Vnames{j,1};
        end
    end
    
    %Permutations Lx and Ly
    %for j=1:nbPerm
    %    p=randperm(size(LxPlot,1));
    %    Lxp=Lx(p,1:size(LxPlot,2));
    %    avgLxp=avg(Lxp);
    %end
end