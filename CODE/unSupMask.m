function mask2=unSupMask(mask,dataClass,numClass)
% function mask2=unSupMask(mask,dataClass,numClass)

%----- 
%----- a recursive function is used
%------
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------   November  2001 ----------------------------
%----------------------------------------------------
%------ input : 
%------         
%------         

%------ no input data is received, error -------------------------
if nargin<1
   help unSupMask;
   mask2=[];
   return;
end;

for counm=1:numClass
    for counc=1:numClass
        pixcoun(counm,counc)=sum(sum((mask==counm)&(dataClass==counc)));
        %pixcoun(counm,counc)=sum(sum((mask==counm)==(dataClass==counc)));
    end
end

%indexCl=1:numClass;
[a,b]=(max(pixcoun));
%indexCl(b(a>0))=[];
[q1,q2]=sort(max(pixcoun));
indexMax=q2(end:-1:1);      %this index points to the Column with higher values
index2=1:numClass;          %this index counts the classes assigned to the new mask
mask2=zeros(size(mask));                 %new mask
indexRem=[];
for counterIndex=indexMax   %the counter is the index to array b
    if ((a(counterIndex>0))&(index2(b(counterIndex))>0))
        %[counterIndex b(counterIndex)]
          
        mask2(mask==b(counterIndex))=counterIndex;
        index2(b(counterIndex))=0;
    else
        indexRem=[indexRem counterIndex];
    end
end

index2(index2==0)=[];
pixcoun=pixcoun(index2,indexRem);
[a,b]=(max(pixcoun));
[q1,q2]=sort(max(pixcoun));
indexMax=q2(end:-1:1);  
for counterIndex=indexMax   %the counter is the index to array b
        mask2(mask==index2(counterIndex))=indexRem(counterIndex);
end

% 
% 
% mask2=mask;
% %b(a==0)=indexCl;
% for counterIndex=indexMax   %the counter is the index to array b
%     if (a(counterIndex>0)&a(b(counterIndex)))
%         index2(counterIndex)=b(counterIndex);
%         index2(b(counterIndex))=counterIndex;
%         a(counterIndex)=0; a(b(counterIndex))=0;
%     end
% end
% %mask3=mask;
% for counm=1:numClass
%     %if a(counm)>0
%     %[ index2(counm) counm]
%         mask2(mask==index2(counm) )=counm;
%         %mask3(mask==counm)=b(counm);
%         %end
% end
% 
% % %[mask2 mask3]
