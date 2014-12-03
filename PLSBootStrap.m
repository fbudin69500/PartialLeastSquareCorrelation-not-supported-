function [Uratio,UConfInf,UConfSup,Vratio,VConfInf,VConfSup]=PLSBootStrap(Data,Groups,splitIndex,U,V,S,sampling,sI)
%sampling with replacement
  nbColumns=size(Data,2);
  Usize=nbColumns-splitIndex;
  Uall=zeros(sampling,Usize,size(sI,2));
  Vall=zeros(sampling,splitIndex,size(sI,2));
  for i=1:sampling
    r=randi([1 size(Data,1)],1,size(Data,1));
    DataSampled=Data(r,1:nbColumns);
    GroupsSampled=Groups(r,1);
    [Xboot,Yboot]=PLSNormalizeData(DataSampled,GroupsSampled,splitIndex);
    [x y]=find(isnan(Yboot));
    if size(x,1) ~= 0
        continue
    end
    [x y]=find(isnan(Xboot));
    if size(x,1) ~= 0
        continue
    end
    Rboot=Yboot'*Xboot;
    %[Uboot,~,Vboot]=svd(Rboot);
    invS=pinv(S);
    Uboot=Rboot*V*invS;
    Vboot=(invS*U'*Rboot)';
    Uall(i,:,:)=Uboot(:,sI);
    Vall(i,:,:)=Vboot(:,sI);
  end
  [Uratio,UConfInf,UConfSup]=PLSConfidenceInterval(U(:,1:size(sI,2)),Uall,sampling,0);
  [Vratio,VConfInf,VConfSup]=PLSConfidenceInterval(V(:,1:size(sI,2)),Vall,sampling,0);
end