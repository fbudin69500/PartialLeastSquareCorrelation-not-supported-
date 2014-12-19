function [NormData]=PLSNormalizeData(data,Groups)
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
end