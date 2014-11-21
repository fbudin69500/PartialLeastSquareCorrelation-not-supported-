function [Mratio,MConfInf,MConfSup]=PLSConfidenceInterval(M,Mall,sampling)
  Mstd=std(Mall);
  Mstdresh=reshape(Mstd,[size(M,1) size(M,2)]);
  Mratio=M./Mstdresh;
  Mmean=mean(Mall);
  Mmeanresh=reshape(Mmean,[size(M,1) size(M,2)]);
  %http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html
  %http://www.losmedanos.edu/Groups/Math/documents/FormulaSheetTable.pdf
  %1.962 corresponds to a confidence interval of 95% with 1000 degrees of
  %freedom
  MConfInf=Mmeanresh-1.962/sqrt(sampling)*Mstdresh;
  MConfSup=Mmeanresh+1.962/sqrt(sampling)*Mstdresh;
end