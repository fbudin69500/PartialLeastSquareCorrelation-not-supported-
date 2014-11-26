function [stableNames]=PLSStable(Mratio,MConfInf,MConfSup,Mnames,index)
supTwo=find(Mratio(:,index)>2);
outsideInf=find( MConfInf(:,index) > 0 );
outsideSup=find( MConfSup(:,index) < 0 );
listStables=[supTwo;outsideInf;outsideSup];
listStablesUnique=unique(listStables);
stableNames=Mnames(listStablesUnique);
end