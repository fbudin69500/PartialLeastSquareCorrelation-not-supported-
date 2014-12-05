function [probaInertia,probaSingularValue,CumulativePercent]=PLSPermutations(X,Y,S,nbPerm,U,V)
    OriginalSingularValues=sum(S,2);
    nbSingularValues=min(size(S));
    OriginalSingularValues(nbSingularValues+1:size(OriginalSingularValues))=[];
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
        [i,s,xp]=MyPermute(X,Y,U,V,S);
        s2=sum(s,2);
        s2(nbSingularValues+1:size(s2))=[];
        InertiaPerm(j)=i;
        SingularValues(:,j)=s2;
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
    SingularValuesVec=reshape(SingularValues,nbSingularValues*nbPerm,1);
    for j=1:nbSingularValues
        LoopSingularValuesVec=SingularValuesVec;
        LoopSingularValuesVec(end+1,1)=OriginalSingularValues(j,1);
        sorted=sort(LoopSingularValuesVec);
        index=find(sorted == OriginalSingularValues(j,1),1 );
        probaSingularValue(1,j)=1-(index-2)/(nbSingularValues*nbPerm);
    end
end