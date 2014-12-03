function [Mratio,MConfInf,MConfSup,MConfInfPerc,MConfSupPerc]=PLSConfidenceInterval(M,Mall,sampling,normal)
  Mstd=std(Mall);
  Mstdresh=reshape(Mstd,[size(M,1) size(M,2)]);
  Mratio=M./Mstdresh;
  Mmean=mean(Mall);
  Mmeanresh=reshape(Mmean,[size(M,1) size(M,2)]);
  %http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html
  %http://www.losmedanos.edu/Groups/Math/documents/FormulaSheetTable.pdf
  %1.962 corresponds to a confidence interval of 95% with 1000 degrees of
  %freedom
  %MConfInf=Mmeanresh-1.962/sqrt(sampling)*Mstdresh;
  %MConfSup=Mmeanresh+1.962/sqrt(sampling)*Mstdresh;
  if normal ~= 0
    MConfInf=Mmeanresh-1.962*Mstdresh;%Mstdresh seems to be the standard error??? so we do not divide by sqrt(N)
    MConfSup=Mmeanresh+1.962*Mstdresh;
  else
    MallSorted=sort(Mall,1);
    MConfInfPerc3D=prctile(MallSorted,2.5,1);
    MConfInf=reshape(MConfInfPerc3D,[size(M,1) size(M,2)]);
    MConfSupPerc3D=prctile(MallSorted,97.5,1);
    MConfSup=reshape(MConfSupPerc3D,[size(M,1) size(M,2)]);
  end
end