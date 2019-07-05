function [featureSpace2,totVal]=sopy3d(data,levelsDown,desiredFeatures)
%function [featureSpace2,totVal]=sopy3d(data,levelsDown,desiredFeatures)
%-----------------------------------------------------------------
%******  Second Orientation PYramid    3D        ******-------------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  5 March 2002    -----------------------------------------
%-----------------------------------------------------------------
%------ input  :  The data  upon which a second orientation  -----
%------           Pyramid will be obtained                   -----
%------           data will be received as spatial domain    -----
%------           optional param: num levels down            -----
%------                           The specific Featues to be
%------                           extracted (can be very slow)
%------ output :  The different levels of the Pyramid        -----
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------





%------ initial revision of input data  --------------------------
if nargin<1   help sopy3d;   featureSpace=[];   return; end;
if nargin==1 levelsDown=1; desiredFeatures=[1:22 24 25 28:32]; end
if nargin==2 
    [r,c,l]=size(levelsDown);
    if (r==1&c==1&l==1) 
        switch levelsDown
        case 1
            desiredFeatures=[1:22 24 25 28:32]; 
        case 2
            desiredFeatures=[1:22 24 25 28:54 56 57 60:64]; 
        case 3
            desiredFeatures=[1:22 24 25 28:54 56 57 60:86 88 89 92:96]; 
        end
    else
        desiredFeatures=levelsDown;
        levelsDown=ceil(max(desiredFeatures)/32);
    end
end

%------ Revise dimensions of data --------------------------------
[rows,cols,levs]=size(data);
%------ !!!!restriction!!! the data sizes must be even.  --------------
if ((rows/4~=floor(rows/4))|(cols/4~=floor(cols/4))|(levs/4~=floor(levs/4)))
    warning('Data dimensions must be divisible by 4');
    featureSpace=[]; help sopy3d; return;
end

%------ Transform the data if needed -----------------------------
if (isreal(data))
    data=fftshift(fftn(data));   
    %data=data;     
end


% tic;
% featureSpace=sopycycle(data,levelsDown,desiredFeatures,rows,cols,levs);  %----- recursive call!!!----
% totalFeatures=size(featureSpace,2);                                     %CELL ARRAY
% totVal=[];
% %------ inverse transform once everything has been set in place -----
% for counterFeatures=1:totalFeatures
%     totVal(counterFeatures)=sum(sum(sum(abs(featureSpace{counterFeatures}))))/rows/cols/levs;
%     if totVal(counterFeatures)~=0
%         tempResult=ifftn(featureSpace{counterFeatures});
%         featureSpace{counterFeatures}=round(abs(tempResult(rows:-1:1,cols:-1:1,levs:-1:1))); %reduces size of feature space
%     end
% end
% if sum(totVal>0)==1      featureSpace=featureSpace{totVal>0};    end
% t1=toc;
% 
% tic;
if length(desiredFeatures)==1
    a=abs(ifftn(filterInThreeD(data,desiredFeatures,rows,cols,levs)));
    if (max(max(max(a)))>100)
        featureSpace2=round(abs(a));
    else
        featureSpace2=(abs(a));           
    end
else    
    for counterFeat=desiredFeatures
        a=abs(ifftn(filterInThreeD(data,counterFeat,rows,cols,levs)));
%        if (max(max(max(a)))>100)
       if (min(min(min(a)))>10)
            featureSpace2{counterFeat}=round(abs(a));
        else
            featureSpace2{counterFeat}=(abs(a));
        end
    end
end
% t2=toc;
% 
% [t1 t2]

%-----------------------------------------------------------------
%------ end of the main function  --------------------------------
%-----------------------------------------------------------------

function dataf=filterInThreeD(dataf,cFeat,rows,cols,levs)

if (cFeat==23|cFeat==26|cFeat==27) dataf=[]; return; end;
if cFeat>32                                     %------ if Feature > 32 then go to next level of pyramid
    cFeat=cFeat-32;                             %------ reposition in the pyramid the feature
    dataf=cDeTessel(filterInThreeD(cTessel(dataf,1/4),cFeat,rows/2,cols/2,levs/2));
    %can be the use od deTessel or clear data and reposition in centre
else
    %------ position of the filter  --------------------------------
    if cFeat==22
        gaussFilt2=gaussF(rows/2,cols/2,levs/2);
        gaussFilt2=gaussFilt2/max(max(max(gaussFilt2 )));
        centralFeature=dataf([1+rows/4:3*rows/4],[1+cols/4:3*cols/4],[1+levs/4:3*levs/4]);           %------ extract region of interest
        dataf=repmat(0,[rows cols levs]);                                                            %------ clears matrix
        dataf([1+rows/4:3*rows/4],[1+cols/4:3*cols/4],[1+levs/4:3*levs/4])=centralFeature.*gaussFilt2;%------ filters and repositions
    else
        levPosition=ceil(cFeat/16);
        colPosition=ceil(cFeat/4)-(levPosition-1)*4;
        rowPosition=mod(cFeat,4);rowPosition(rowPosition==0)=4;
        gaussFilt1=gaussF(rows/4,cols/4,levs/4);
        gaussFilt1=gaussFilt1/max(max(max(gaussFilt1 )));
        rIndex=  positionIndex(rows,rowPosition);
        cIndex=  positionIndex(cols,colPosition);
        lIndex=  positionIndex(levs,levPosition);
        %----- For 3D it is important to consider both sides, since they can be a different levels    
%         rIndex2=  positionIndexM(rows,rowPosition);
%         cIndex2=  positionIndexM(cols,colPosition);
%         lIndex2=  positionIndexM(levs,levPosition);
        centralFeature=dataf([rIndex],[cIndex],[lIndex]);
%         centralFeatureM=dataf([rIndex2],[cIndex2],[lIndex2]);     %mirror image
        dataf=repmat(0,[rows cols levs]);                                    %------ clears matrix
        dataf([rIndex],[cIndex],[lIndex])=centralFeature.*gaussFilt1;        %------ filters and repositions
%         dataf([rIndex2],[cIndex2],[lIndex2])=centralFeatureM.*gaussFilt1;    %------ filters and repositions the mirror
    end
end
% % %-----------------------------------------------------------------

function rowIndex=positionIndex(rows,rPosition)
switch rPosition
case 1
    rowIndex=1+3*rows/4:rows;
case 2
    rowIndex=1+rows/2:3*rows/4;
case 3
    rowIndex=1+rows/4:rows/2;
case 4
    rowIndex=1:rows/4;          %--- rows ----
end    
% % %-----------------------------------------------------------------
function rowIndex=positionIndexM(rows,rPosition)      %Will mirror the previous function
switch rPosition
case 1
    rowIndex=1:rows/4;
case 2
    rowIndex=1+rows/4:rows/2;
case 3
    rowIndex=1+rows/2:3*rows/4;
case 4
    rowIndex=1+3*rows/4:rows;          %--- rows ----
end    




