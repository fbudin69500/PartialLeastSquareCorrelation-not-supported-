function Mnorm=MyNormalizedMatrix(M)
    averageVal=mean(M,1);
    averageMat=ones(size(M,1),1)*averageVal;
    %stdVal=std(M);
    %stdMat=ones(size(M,1),1)*stdVal;
    %Mnorm=(M-averageMat)./stdMat;
    MnormAv=M-averageMat;
    Mnorm=normc(MnormAv);
end