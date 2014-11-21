function [X,Y]=PLSNormalizeData(data,Groups,splitIndex)
    rowsSize=size(data,2);
    nbgroups=max(Groups);
    for i=1:nbgroups
        index=find( Groups == i);
        M=data(index,:);
        NormGroup=MyNormalizedMatrix(M);
        if exist('NormData','var')
            NormData=[NormData;NormGroup];
        else
            NormData=NormGroup;
        end
    end
    X=NormData(:,1:splitIndex);
    Y=NormData(:,splitIndex+1:end);
end