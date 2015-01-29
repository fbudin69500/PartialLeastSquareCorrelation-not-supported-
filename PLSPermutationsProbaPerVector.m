function [probaInertia,probaSingularValue,CumulativePercent]=PLSPermutationsProbaPerVector(X,Y,S,nbPerm,Groups)
    OriginalSingularValues=sum(S,2);%Transform square matrix to vector
    nbSingularValues=min(size(S));%compute vector length
    OriginalSingularValues(nbSingularValues+1:size(OriginalSingularValues))=[];%Resize Singular value vector to remove trailing 0's
    TotSq=sum(OriginalSingularValues.^2);
    %percentage of variance explained by the correspond-ing pairs of latent variables
    Percent=OriginalSingularValues.^2./TotSq;
    CumulativePercent=cumsum(Percent);
    %Inertia
    Inertia=sum(OriginalSingularValues);
    %Permutations
    InertiaPerm=zeros(1,nbPerm);
    SingularValues=zeros(nbSingularValues,nbPerm);
    for j=1:nbPerm
        [i,s,xp]=MyPermute(X,Y,Groups);
        if i < -1
            j=j-1;
        else
          s2=sum(s,2);
          s2(nbSingularValues+1:size(s2))=[];
          InertiaPerm(j)=i;
          SingularValues(:,j)=s2;
        end
    end
    %Inertia of permutations
    %hist(InertiaPerm)
    %percentile=prctile(InertiaPerm,95)
    InertiaPerm(end+1)=Inertia;
    sorted=sort(InertiaPerm);
    index=find(sorted == Inertia );
    probaInertia=1-(index-2)/nbPerm;
    %Significance of singular values
    probaSingularValue=zeros(1,nbSingularValues);
    for j=1:nbSingularValues
        %LoopSingularValuesVec=SingularValuesVec;
        LoopSingularValuesVec=SingularValues(j,:);
        LoopSingularValuesVec(1,end+1)=OriginalSingularValues(j,1);
        sorted=sort(LoopSingularValuesVec);
        index=find(sorted == OriginalSingularValues(j,1),1 );
        probaSingularValue(1,j)=1-(index-2)/nbPerm;
    end
end