function [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,Groups,splitIndex,U,V,S,sampling,sI)
%sampling with replacement
  nbColumns=size(Data,2);
  Uall=zeros(sampling,size(U,1),size(sI,2));
  Vall=zeros(sampling,size(V,1),size(sI,2));
  M=U*S*V;
  for i=1:sampling
    r=randi([1 size(Data,1)],1,size(Data,1));
    DataSampled=Data(r,1:nbColumns);
    GroupsSampled=Groups(r,1);
    normedData=PLSNormalizeData(DataSampled,GroupsSampled);
    %Split data
    Xboot=normedData(:,1:splitIndex);
    Yboot=normedData(:,splitIndex+1:end);
    RbootT=StackGroups(Xboot,Yboot,Groups);
    [x y]=find(isnan(Yboot));
    if size(x,1) ~= 0
        continue
    end
    [x y]=find(isnan(Xboot));
    if size(x,1) ~= 0
        continue
    end
    %RbootT=Yboot'*Xboot;
    %rotation correction
    [U2,S2,V2]=svd(RbootT'*M);
    R=U2*V2';
    Rboot=RbootT*R;
    invS=pinv(S);
    Uboot=Rboot*V*invS;
    Vboot=(invS*U'*Rboot)';
    Uall(i,:,:)=Uboot(:,sI);
    Vall(i,:,:)=Vboot(:,sI);
  end
  [Uratio,UConfInf,UConfSup]=PLSConfidenceInterval(U(:,1:size(sI,2)),Uall,sampling,0);
  [Vratio,VConfInf,VConfSup]=PLSConfidenceInterval(V(:,1:size(sI,2)),Vall,sampling,0);
end