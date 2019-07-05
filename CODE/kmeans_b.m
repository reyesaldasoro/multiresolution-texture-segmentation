function [dataClass,err,means,err2]=kmeans_b(means,mask,numIterations,varargin)  
%function [dataClass,err,means]=kmeans_b(means,mask,numIterations,varargin)  
%-----------------------------------------------------------------
%******   K-means supervised/unsupervised clustering      ******-------------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  18 September  2003 ---------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------
%-----------------------------------------------------------------
%------ input  :  means / MASKS / ITERATiONS / DATA                       -----  
%------           means can be a value (numClasses) or the prototypes
%------           means can be supplied (from LBG for instance) or generated randomly
%------           Data and Mask can be 2d/3d mask is always required
%------ output :  The data classified into the classes predefined   
%------           The error against the mask
%-----------------------------------------------------------------
%------ The parameters of the algorithm are:                  ----
%------ dataSets        : the number of sets as input (read from data)---
%------ numIterations   : number of iterations of the process    ---
%------ numClasses      : number of Classes for classification    ----
%------ means           : the means of the predefined classes     ----
%----- SEE ALSO   minDist.m  arrangeData.m
%-----------------------------------------------------------------
%------ no input data is received, error -------------------------
if nargin<4   help kmeans;    dataClass=[];  err=[]; return; end;
if numIterations<1 numIterations=1; end
if ~exist('mask') error('mask is not provided! It is needed to determine size of data output (could be only ones)'); end
% %------------------------------------------------------------------------
% %------ rearrange data, everything into column form  for 2D /3D cases ---
% %------------------------------------------------------------------------
[data,rows,cols,levs,numFeats,RCL]=arrangeData(varargin);
%------ stdData will keep standard deviation of the data per level
%------ it is used to normalise the data and the distances
deltaV=1e-10; %alpha2=0.8;
meanData=mean(data);
stdData=deltaV+std(data);
%----- mask ------
totmask=sum(((mask(:)>0)));   
[rowMask,colMask,levMask]=size(mask);
RCLM=rowMask*colMask*levMask;        %----- size of data sets rearranged
if RCLM~=RCL  warning('Mask and Data Dimensions are not consistent');    dataClass=[]; err=[]; return; end
mask=reshape(mask,RCL,1);

%------ means ------
if all(size(means)==1)
    %----- warning('means not present, a non-supervised approach will follow!'); 
    %numFeats=prod( size(data))/prod(size(mask));
    numClasses=means;
    %----- generate random means- get the characteristics of the data in each feature (dimension)
    minData=min(data);     maxData=max(data);     
    range=maxData-minData;
    %%------ random (uniform) points within  [min,max] of each dimension
    %means=(range'*ones(1,numClasses)).*rand(numFeats,numClasses)+(minData'*ones(1,numClasses))
    %------ random (normal) points within  [miu,std] of each dimension        in general randn is better
    %means=(meanData'*ones(1,numClasses))+(stdData'*ones(1,numClasses)).*randn(numFeats,numClasses);
    %------ use LBG algorithm (good but slower)
    means=LBG(data,numClasses);   %means=LBG(data',numClasses);
    [numFeatsM,numClasses]=size(means);
    %numFeatsM=numFeats;
else
    [numFeatsM,numClasses]=size(means);
end
% % %------ last check of dimensions
if (numFeatsM~=numFeats)
     warning('Number of Features is not consistent with dimensions of means');    dataClass=[]; err=[]; return; 
end

%------ input process is finished, the variables as column vectors are:
%------         data     mask          means

%------ random classes (uniform) for every point
dataClass= 1+(numClasses-1)*rand(RCL,1);
%dataClass=zeros(RCL,1);


%-----------------------------------------------------------------
%------ The k-means cycle begins here ----------------------------
%-----------------------------------------------------------------
for counter=1:numIterations
    err=[];
    %----- determine minimum distance to means/prototypes and assign to class
    dataClass2=minDist(data,means);
    %-----  measure error if mask is provided
    if exist('mask')    
        err=[err;sum(dataClass2==mask)/totmask];  
        %totBone=sum(sum(sum(mask>0)));totNoBone=sum(sum(sum(mask==0)));
        %dd=dataClass2;boneC=2;r2=strcat('### ',num2str([sum(sum(sum((mask>0)&(dd==boneC))))/totBone sum(sum(sum((mask==0)&(dd==boneC))))/totNoBone sum(sum(sum((mask>0)&(dd~=boneC))))/totBone sum(sum(sum((mask==0)&(dd~=boneC))))/totNoBone]))
    end
    %------ determine the number of pixels that have changed position since the last classification
    pixChanged=sum(dataClass2~=dataClass);
    dataClass=dataClass2;   %used to keep track of changes
   %------ exit if the number of pixels that change is less than 0.1% of total
    if pixChanged<(0.001*rows*cols) 
        dataClass=reshape(dataClass2,rowMask,colMask,levMask); return;
    end;
    %------ recalculate class will distribute dataClass in several slices of dataC
    for counter2=1:numClasses
        inClass=(dataClass==counter2);
        totClass=sum(inClass);
        if totClass>0
            dataClassX=data(inClass,:);
            %new means can be weigthed by the number of pixels added to the new class
            means(:,counter2)=((RCL-totClass)*means(:,counter2)+(totClass)*mean(dataClassX)')/RCL;
            % means(:,counter2)=mean(dataClassX)';
        else
            %----- we are in a DEAD Class, no pixels assigned to it, what to do????
            %----- try reassigning a new value at random
            means(:,counter2)=meanData'+(stdData').*randn(numFeats,1);
        end
    end
    %------ contiguity and compactness
  %   CC=contiguity(reshape(dataClass2,rows,cols,levs),numClasses);
  %    [rc,tc,cv]=compactness(data,dataClass2,numClasses);
    %------- just print the stage process
   % [counter err pixChanged/rows/cols CC mean2(rc) ]
    %figure;
    %surfdat(reshape(dataClass2,rowMask,colMask,levMask));
end     
%
%-----------------------------------------------------------------
%------ End of the k-means cycle ---------------------------------
%-----------------------------------------------------------------
dataClass=reshape(dataClass2,rowMask,colMask,levMask);

if nargout>3
    mask=reshape(mask,rowMask,colMask,levMask);

    %%%% addition 7 may 2004 in unsupervised cases, mask is not reliable  so,
    %%%% generate that mask to minimise the misclassification
    mask2=unSupMask(mask,dataClass,numClasses);
    err2=[sum(sum(sum(dataClass==mask2)))/totmask];
end