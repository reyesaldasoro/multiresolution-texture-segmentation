function  [dataClass,errClass]=multiLClass(data,meansC,mask,levsP,algorMethod,dampP)
%function  [dataClass,errClass]=mVts(data,meansC,mask,levsP,algorMethod,dampP)
%---------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro   --------------------
%------             PHD     the University of Warwick   --------------------
%------  Supervisor :   Abhir Bhalerao    ----------------------------------
%------  8 September 2003 --------------------------------------------------
%---------------------------------------------------------------------------
%------ For a description and explanation please refer to:------------------
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts -------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%------ Input : data -      the data (measurement space) to be classified, 
%------                     3D in case of an image, 4D in case of a volume
%------         meansC      The means of the classes to be segmented OR the number of classes
%------         mask        The mask to measure the error (correct classification actually)
%------         levsP       number of levels of a Pyramid to climb
%------         algorMethod Option for the butterfly filters, default 4 mode-mode
%------         dampP       Damping parameter, how close to mean in butterfly filter
%------                     Refer to figure 5.8 of the thesis


 
global dampParam numFeats numClasses butterf levsPyr means rowsM colsM levsM algor alpha1 beta1

%----- regular dimensions check ----------------------
%----- original data to classify can be 2D / 3D ===>
%----- which implies data (Feature Space) is 3D /4D --
[rows,cols,levs,dataSets]=size(data); 
[rowsM,colsM,levsM]=size(mask);
if ((rows~=rowsM)|(cols~=colsM))    warning('mask and data dimensions are not consistent');     dataout=[]; return;end
if xor(dataSets==1,levsM==1)        warning('mask and data dimensions are not consistent 2D/3D');     dataout=[]; return;end
if ~exist('meansC')' 
    warning('Either Means or Number of Classes is required');     dataout=[]; return;
else
    if all(size(meansC)==1)
        %warning('means not present, a non-supervised approach will follow!'); 
        numFeats=prod( size(data))/prod(size(mask));
        numClasses=meansC;
        %means=ones(numFeats,numClasses);%
        %means=LBG(data,numClasses);
        means=meansC;
    else
        means=meansC; 
        [numFeats,numClasses]=size(means);
    end
end


%------- parameters of the algorithm
alpha1=1;beta1=2;               %------- found to optimise misclassification for Randen image   %alpha1=alpha2;beta1=beta2;
if ~exist('levsP')          levsPyr=log2(rows)-1;       else    levsPyr=levsP;      end
%---- this dampens the new position, factor of the distance between classes
if ~exist('dampP')          dampParam=-1.0;             else    dampParam=dampP;    end
if ~exist('algorMethod')    algor=4;                    else    algor=algorMethod;  end


%----- create butterflies 5 x 5 as default  ----------------------
if levsM==1;                levsB=1;                    else    levsB=5;            end
butterf=createButter(5,5,levsB,4);

%------ if filter the data in the SOP is needed
%levsSOP=floor(numFeats/7);
%if levs==1
%    dataSOP=sopy(data,levsSOP,1);    %----- use 0 for complex result 1 for mag
%elseif levs==numFeats
dataSOP=data;
%end

%----- call to the recursive process    
[dataClass,boundary,orientationFi,orientationTheta,errClass]=classifPyr (dataSOP,mask,1);
clear  global dampParam numFeats numClasses butterf levsPyr means rowsM colsM levsM algor
%----------------------------------------------------
%----- End of main routine ----------------------------
%----------------------------------------------------

%----------------------------------------------------
%---------- REcursive call process      ------------------
%----------------------------------------------------
function [dataClass,boundary,orientationFi,orientationTheta,errClass]= classifPyr(qtdata,maskPyr,presentLevel);
global dampParam numFeats numClasses butterf levsPyr means rowsM colsM levsM
[rPyr,cPyr,lPyr,numFeats2]=size(qtdata);
%------------ if top level reached, classify and propagate downwards
if levsPyr==presentLevel
    %------- LBG linde-buzo-gray initial values:
    if all(size(means)==1)
        means=LBG(qtdata,numClasses);
    end
     [dataArr,rows2,cols2,levs2,numFeats2,RCL2]=arrangeData(qtdata);
    %-------- following lines are for the spatial references xx yy zz
    maxData=max(max(dataArr)); minData=min(min(dataArr));
    xx= repmat(ones(rows2,1)*[0:cols2-1],[1 1 levs2]);
    yy= repmat([0:rows2-1]'*ones(1,cols2),[1 1 levs2]);
    for counterZ=1:levs2 zz(:,:,counterZ)=ones(rows2,cols2)*counterZ; end
    
    %   xx=xx(:);yy=yy(:);zz=zz(:);             %only turn into column if hierarchical is used
    xxn=minData+(maxData-minData)*[xx/cols2];
    yyn=minData+(maxData-minData)*[yy/rows2];
    zzn=minData+(maxData-minData)*[zz/rows2];
    %xx=ones(rows2,1)*[1:cols2];yy=[1:rows2]'*ones(1,cols2);xx=xx(:);yy=yy(:);
    %xx=0.01*randn(rows2,cols2)+ ones(rows2,1)*[0:cols2-1];yy=0.01*randn(rows2,cols2)+[0:rows2-1]'*ones(1,cols2);xx=xx(:);yy=yy(:);
    
    %------- next 5 lines are used in case agglomerative clustering is required  (good results but requires a lot of memory
    %------= for the distance!
    %dataArr=normData(dataArr);
    %------- next two lines normalise the data so that every feature has STD=1 ----
   
    %------ using MST to cluster  there is something wrong .... the dataClass looks awful :-( check later
            %[q1,q2,distanceData]=MST([ dataArr xxn yyn zzn ]);
            %[q1,q2,distanceData]=MST([dataArr]);
            %classDat=clusterMST(q1,q2,numClasses);
            %     dataClass{presentLevel}=reshape(classDat,rows2,cols2,levs2);
            %     for counClasses=1:numClasses   means2(:,counClasses)= mean(dataArr(classDat==counClasses,:))'; end
             %------next lines are using aggclust procedure, does not work for higher dimensions, MEMORY problems
% %     [level,classDat]=aggclust(distanceData);
% %     dataClass{presentLevel}=reshape(classDat(:,RCL2-numClasses+1),rows2,cols2,levs2);
% %      for counClasses=1:numClasses  
% %         means2(:,counClasses)= mean(dataArr(level(RCL2-numClasses+1).cluster{counClasses},:))'; 
% %         %sigmas2(:,counClasses)= std(dataArr(level(RCL2-numClasses+1).cluster{counClasses},:))'; 
% %     end
% %     means=means2;
    
    %------- next line is used for kmeans clustering
    [dataClass{presentLevel},err,means]=kmeans_b(means,maskPyr,1,qtdata);errClass=err(end);         % this does not includes xx,yy,zz
    %qtdata(:,:,:,numFeats2+1)=xx;qtdata(:,:,:,numFeats2+2)=yy;qtdata(:,:,:,numFeats2+3)=zz;        % append xx yy zz for 3d data
    %[dataClass{presentLevel},err,means]=kmeans_b(numClasses,maskPyr,1,qtdata);errClass=err(end);   % this does not includes xx,yy,zz
    %qtdata=qtdata(:,:,:,[1:numFeats2]);
    means=means(1:numFeats2,:);
    %figure(2);surfdat(dataClass{presentLevel});%beep    %err
     %------- next lines are used to get new means
   
    
   [dataClass{presentLevel},boundary,orientationFi,orientationTheta]=refineData(dataClass{presentLevel},numClasses);%,rPyr,cPyr,4);
    %figure(2);surfdat(dataClass{presentLevel});
    errClass=sum(sum(sum(dataClass{presentLevel}==maskPyr)))/prod(size(maskPyr)); 
else
    %------------ if top level NOT reached, reduce and push upwards---------------------------
    qtdata2=upQuadTree(qtdata);
    [dataClass,boundary,orientationFi,orientationTheta,errClass]= classifPyr(qtdata2,maskPyr(1:2:end,1:2:end,1:2:end),presentLevel+1);
    %-----  classification of the highest level is stored in *dataClass* ---------------------------
    %--------------------propagation to lower levels-----------------------------------------
    dataClass{presentLevel}=expandu(dataClass{presentLevel+1});     
    boundary=expandu(boundary);       
    orientationFi=expandu(orientationFi);
    orientationTheta=expandu(orientationTheta);
    %-------------------- markov random fields -----------------------------------------
    %dataClass{presentLevel} =  markovFields(qtdata,dataClass{presentLevel},7);
    %--------------------butterfly filters and refinement-----------------------------------------
    if prod(size(dataClass{presentLevel}))> (64*64*18)
        dataClass{presentLevel} =butterFilter2(qtdata,dataClass{presentLevel},boundary,orientationFi,orientationTheta,rPyr,cPyr,lPyr);
        [dataClass{presentLevel},boundary,orientationFi,orientationTheta]=refineData(dataClass{presentLevel},numClasses);%,rPyr,cPyr,4);
    else
        dataClass{presentLevel} =butterFilter(qtdata,dataClass{presentLevel},boundary,orientationFi,orientationTheta,rPyr,cPyr,lPyr);
        [dataClass{presentLevel},boundary,orientationFi,orientationTheta]=refineData(dataClass{presentLevel},numClasses);%,rPyr,cPyr,4);
    end
    %--------------------   calculate errors        -----------------------------------------
    errClass=sum(sum(sum(dataClass{presentLevel}==maskPyr)))/prod(size(maskPyr))   ;
    %errClass2=sum(sum(sum(dataClass2{presentLevel}==maskPyr)))/prod(size(maskPyr))   ;
    [errClass ]
    subplot(312);surfdat(dataClass{presentLevel});drawnow;
    %subplot(313);surfdat(dataClass2{presentLevel})
    %dataClass{presentLevel}=dataClass2{presentLevel};
end
%----------------------------------------------------
%------------ End of cycle----------------------------
%----------------------------------------------------



%----------------------------------------------------
%------------ subfunction  markovFields  --------------
%------------ reclassify depending on the neighbours ---
%----------------------------------------------------
function  dataClass =  markovFields(qtdata,dataClass,betaVal);

if ~exist('betaVal') betaVal=1; end
global  numClasses 
 [dataArr,rows2,cols2,levs2,numFeats2,RCL2]=arrangeData(qtdata);
 dClass=reshape(dataClass,RCL2,1);
%------------------get means and stds ------------------------
 sigmas=zeros(numFeats2,numFeats2,numClasses);
for cClass=1:numClasses  
        tempDat=(dataArr(dClass==cClass,:));
        means(:,cClass)= (mean(tempDat))';
        sigmas(:,:,cClass)= diag(std(tempDat));
end
%------------------get number of neighbour elements in class ---
dClass2=vol2col(padData(dataClass,1),[3 3 3]);
[numNeighbs,numElements]=size(dClass2);
dClass2(1+floor(numNeighbs/2),:)=[];

%------------------calculate normal probability and potentials ------------------------
for cClass=1:numClasses 
    Vc=betaVal*(sum(dClass2==cClass))';
    aPrioriP(:,cClass)=exp(-(numNeighbs-1)- (2*Vc));
    for cElem=1:1024:numElements
        conditionalProb(cElem:cElem+1023,cClass)=normProb(dataArr(cElem:cElem+1023,:),means(:,cClass),sigmas(:,:,cClass));
    end
end

[dataMinim,newClass]=min(aPrioriP.*conditionalProb,[],2);
dataClass=reshape(newClass,rows2,cols2,levs2);



%----------------------------------------------------
%------------ subfunction  upQuad Tree  ------------
%------------ climb one level in the quad tree 2D/3D
%--------------------------------------------------- 
function qtdata2=upQuadTree(qtdata);
global numFeats
[rQT,cQT,lQT,fQT]=size(qtdata);
if size(qtdata,4)==1        %---- 2D Data
    qtdata2=zeros(rQT/2,cQT/2,lQT);
    for cSOP=1:numFeats      qtdata2(:,:,cSOP)=reduceu(qtdata(:,:,cSOP));       end
else
    qtdata2=zeros(rQT/2,cQT/2,lQT/2,fQT);
    for cSOP=1:numFeats      qtdata2(:,:,:,cSOP)=reduceu(qtdata(:,:,:,cSOP));   end
end

%----------------------------------------------------
%------------ subfunction  butterFilter  ------------
%------------ filter the relevant data witha butterfly filter
%--------------------------------------------------- 
function    [dataClass] =  butterFilter (qtdata,dataClass,boundary,orientationFi,orientationTheta,rPyr,cPyr,lPyr)
%------ input :     qtdata data to be classified
%------             means - prototype of the classes
%------             dataClass las classification
%------             boundary - actual pixels to be filtered 
%------ output:     smooth classification
%------ implicit variables : butterf - the filter to be used for every boundary pixesl

%----- first check if there IS a boundary, if not there is only one class
%----- and there is no need of refinement
if sum(sum(sum(boundary)))>0 
    delta=1e-10;
    global dampParam numFeats numClasses butterf levsPyr means algor alpha1 beta1
    %----- dimensions of the variables ------------------------------------------- 
    isData3D=size(dataClass,3);             sizButterf=size(butterf);               
    xSizeButterf=sizButterf(1);             numNeighbours=prod(sizButterf(1:3));
    boundaryCol=boundary(:);                boundaryPixels=sum(boundaryCol~=0);
    numPadPixels=(xSizeButterf-1)/2;

    %----- if a point is not part of the boundary, take its class mean value (prototype)
    if isData3D>1        reppattern=[1 1 1 numFeats];      reppattern2=[rPyr cPyr lPyr 1];
    else        reppattern=[1 1 numFeats];        reppattern2=[rPyr cPyr  1];
    end
    boundaryMultiD=repmat(boundary,reppattern);
    dataClassMultiD=repmat(dataClass,reppattern);
    qtdata(boundaryMultiD==0)=0;
    for cClass=1:numClasses
        qtdata=qtdata+(boundaryMultiD==0).*(dataClassMultiD==cClass).*repmat(reshape(means(:,cClass),reppattern),reppattern2);
    end
    
        %----- pad  with SAME values on the edges
        %----- size of butterfly 3, 5, 7 2n+1 =>  1, 2, 3  ...n pad values
        dataClass3=padData(dataClass,numPadPixels);
        %qtdata3=padData(qtdata,numPadPixels);
        qtdata3=padData(qtdata,numPadPixels,size(dataClass3));
        %---- convert into rows
        dataClassCol=vol2col(dataClass3,xSizeButterf*[1 1 1],boundary);
        if  isData3D>1     
            for cFeats=1:numFeats qtdataCol(:,:,cFeats)=vol2col(qtdata3(:,:,:,cFeats),xSizeButterf*[1 1 1],boundary); end
            orientationTheta=orientationTheta(:);
            orientationTheta(boundaryCol==0)=[];
        else
            for cFeats=1:numFeats qtdataCol(:,:,cFeats)=im2colRed(qtdata3(:,:,cFeats),xSizeButterf*[1 1],boundary); end; 
        end
        %----- now reduce every thing to the boundary pixels
        %dataClassRefined=dataClassCol(ceil(numNeighbours/2),:);
        %qtdataCol(:,boundaryCol==0,:)=[];%dataClassCol(:,boundaryCol==0,:)=[];%
        dataClassRefined=dataClass(:);
        orientationFi=orientationFi(:);         orientationFi(boundaryCol==0)=[];
        %----- Construct a butterfly filter mask depending on the orientation
        butterMask=zeros(size(dataClassCol));
        if isData3D>1
            for cFi=1:4
                indexFi=(orientationFi==cFi);
                for cTheta=1:4
                    indexTheta=(orientationTheta==cTheta);
                    butterMask(:,indexFi&indexTheta)=reshape(butterf(:,:,:,cFi,cTheta),numNeighbours,1)*ones(1,sum((orientationFi==cFi)&(orientationTheta==cTheta)));
                end
            end
        else
            for cOrient=1:4
                butterMask(:,orientationFi==cOrient)=reshape(butterf(:,:,:,cOrient),numNeighbours,1)*ones(1,sum(orientationFi==cOrient));
            end
        end
        
        bMask1=(butterMask==1);
        bMask3=(butterMask==3);
        
        %----- actual filtering process
        [q1,mod1]=mode(dataClassCol.*(bMask1),2);
        [q3,mod3]=mode(dataClassCol.*(bMask3),2);
        dClass1=(dataClassCol==repmat(mod1,numNeighbours,1));
        dClass3=(dataClassCol==repmat(mod3,numNeighbours,1));
        %------ UN-normalised Distance to the means of the modes Class --------------------------
        
        %---- new variable algor to determine which options to use all / mode
        %----  1 -   all / all
        %----  2 -   all / mode
        %----  3 -   mode / all
        %----  4 -   mode / mode
        
        if (algor==1)|(algor==2)
            %------ new points are calculated  from all points in the butterfly
            numPixAver1=sum(bMask1);
            numPixAver3=sum(bMask3);
            miu1=sum(qtdataCol.*repmat(bMask1,[1 1 numFeats]))./repmat(numPixAver1,[1 1 numFeats]);
            miu3=sum(qtdataCol.*repmat(bMask3,[1 1 numFeats]))./repmat(numPixAver3,[1 1 numFeats]);
            %-------- the following lines are used in case a per pixel dampParam is
            %-------- required, otherwise keep a fixed value for the whole distribution
            miusq1=sum((qtdataCol.^2).*repmat(bMask1,[1 1 numFeats]))./repmat(numPixAver1,[1 1 numFeats]);
            miusq3=sum((qtdataCol.^2).*repmat(bMask3,[1 1 numFeats]))./repmat(numPixAver3,[1 1 numFeats]);
            var1=miusq1-miu1.^2;
            var3=miusq3-miu3.^2;
            damp0=((miu1-miu3).^2)./(var1+var3);
            dampParam=1./(1+exp(1*(5-1*damp0)));
            r1=  dampParam.*qtdataCol(ceil(numNeighbours/2),:,:) +  (1-dampParam).*miu1;
            r3=  dampParam.*qtdataCol(ceil(numNeighbours/2),:,:) +  (1-dampParam).*miu3;
        elseif (algor==3)|(algor==4)
            %------ new points are calculated only from those points belonging to the modes
            numPixAver1=sum(dClass1&bMask1);
            numPixAver3=sum(dClass3&bMask3);
            %-------- the following lines are used in case a per pixel dampParam is
            %-------- required, otherwise keep a fixed value for the whole distribution
            miu1=sum(qtdataCol.*repmat(dClass1&bMask1,[1 1 numFeats]))./repmat(numPixAver1,[1 1 numFeats]);
            miu3=sum(qtdataCol.*repmat(dClass3&bMask3,[1 1 numFeats]))./repmat(numPixAver3,[1 1 numFeats]);
            miusq1=sum((qtdataCol.^2).*repmat(bMask1,[1 1 numFeats]))./repmat(numPixAver1,[1 1 numFeats]);
            miusq3=sum((qtdataCol.^2).*repmat(bMask3,[1 1 numFeats]))./repmat(numPixAver3,[1 1 numFeats]);
            var1=miusq1-miu1.^2;
            var3=miusq3-miu3.^2;
            damp0=((miu1-miu3).^2)./(var1+var3+delta);
            if (any(size(dampParam)~=1))|(dampParam<0)
                %             dampParam=1./(1+exp(1*(5-1*damp0)));
                dampParam=1./(1+exp(alpha1*(beta1-1*damp0)));
            end
            
            r1=dampParam.*qtdataCol(ceil(numNeighbours/2),:,:) + (1-dampParam).*sum(qtdataCol.*repmat(dClass1&bMask1,[1 1 numFeats]))./repmat(numPixAver1,[1 1 numFeats]);
            r3=dampParam.*qtdataCol(ceil(numNeighbours/2),:,:) + (1-dampParam).*sum(qtdataCol.*repmat(dClass3&bMask3,[1 1 numFeats]))./repmat(numPixAver3,[1 1 numFeats]);
        end
        %----- the squeeze reshapes the new points in the format
        %-----   [numBoundaryPoints numFeats]
        r11=squeeze(r1); r33=squeeze(r3);
        if (algor==1)|(algor==3)
            %----compare between all possible classes, 
            varPoints1=repmat(var(r11),[boundaryPixels 1 numClasses]);
            varPoints3=repmat(var(r33),[boundaryPixels 1 numClasses]);
            meansExp=repmat(reshape(means,[1 numFeats numClasses]),[ boundaryPixels 1 1]);
            dist1=sum(((repmat(r11,[1 1 numClasses])-meansExp).^2)./varPoints1,2);
            dist3=sum(((repmat(r33,[1 1 numClasses])-meansExp).^2)./varPoints1,2);
            [d1min,pixClass]=min(min([dist1 dist3],[],2),[],3);
        elseif (algor==2)|(algor==4)
            %----compare ONLY between two classes, those which  are most common on each butterfly
            %for normalised distances divide  by the variance of the POINTS 
            varPoints1=repmat(var(r11),[boundaryPixels 1 1]);
            varPoints3=repmat(var(r33),[boundaryPixels 1 1]);
            dist11=sum(((r11-means(:,mod1)').^2)./varPoints1,2);
            dist13=sum(((r11-means(:,mod3)').^2)./varPoints1,2);
            dist31=sum(((r33-means(:,mod1)').^2)./varPoints3,2);
            dist33=sum(((r33-means(:,mod3)').^2)./varPoints3,2);
            %%%----- un-normalised distance ---------------
            %      dist11=sum(((r11-means(:,mod1)').^2),2);
            %      dist13=sum(((r11-means(:,mod3)').^2),2);
            %      dist31=sum(((r33-means(:,mod1)').^2),2);
            %      dist33=sum(((r33-means(:,mod3)').^2),2);
            [dataMinim1,pixClass1]=min([dist11';dist13';dist31';dist33']); 
            %------ class assignment --------------------
            pixClass=pixClass1;
            pixClass(pixClass1==1)=mod1(pixClass1==1);
            pixClass(pixClass1==3)=mod1(pixClass1==3);
            pixClass(pixClass1==2)=mod3(pixClass1==2);
            pixClass(pixClass1==4)=mod3(pixClass1==4);
            
        end
        dataClassRefined(boundaryCol~=0)=pixClass;  %!!!! it was pixClass1 !!!!
        dataClass=reshape(dataClassRefined,size(dataClass));

end
%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------
%------------ subfunction refineData  -----------------------------------------------------------
%------------ compares against neighbours and determines if it is a booundary or not and --------
%------------ then erases the isolates pixels and obtains the orientations-----------------------
%------------------------------------------------------------------------------------------------

function [dataClass,boundary,orientationFi,orientationTheta]=refineData(dataClass,numClasses)
[rPyr,cPyr,lPyr]=size(dataClass);
numOrient=4;
[dClassCompared,neighVal]=pixvsn(dataClass,numClasses);             %----- this pixel comparison will only change isolated pixels
%[dataClass,neighVal]=pixvsn(dataClass,[3 3],numClasses);           %----- this pixel comparison will change pixels according to the rules in pixvsn
%------ to change pixels whose 8 neighbours are different from themselves
dataClass(neighVal==1)=(dClassCompared(neighVal==1));
%-----  or alternatively change those whose 7,6 ... neighbours are different
%dataClass(neighVal<3)=(dClassCompared(neighVal<3));

%----- boundaries are established from the REASSIGNED matrix, therefore
%----- the neighbours values is calculated again
[dClassCompared,neighVal]=pixvsn(dataClass,numClasses);       %----- this pixel comparison will only change isolated pixels
%---- the boundary is established by all those pixels with neighbours different from themselves
[boundC,bound]=pixvsn(ones(size(dataClass)),1);
boundary=(neighVal<bound);
orientationFi=zeros(size(boundary));

%----- orientation values ---------------------------------------------------------
difRows=sign(diff(dataClass,1,1));                      % |   vertical difference    the sign is the only
difCols=sign(diff(dataClass,1,2));                      % -   horizontal difference  thing important
or3=([zeros(1,cPyr,lPyr);-difRows]+[difRows;zeros(1,cPyr,lPyr)] +i*([zeros(rPyr,1,lPyr) -difCols]+[difCols zeros(rPyr,1,lPyr)]));
%----- if the data is big enough then a blur can help the classification
if size(dataClass,1)>16 
    sizeFiltG=3;                                    % just spread a little
    filtg=gaussF(sizeFiltG,sizeFiltG,1);
    if lPyr>1
        for counterL=1:lPyr   or3(:,:,counterL)=conv2(or3(:,:,counterL),filtg,'same');        end
    else
        or3=conv2(or3,filtg,'same');                    
    end
end
%----- orientation matrix obtained
or4=(angle(or3)+pi*(angle(or3)<0))/pi;              % adding pi corrects the negative values   0<or4<1
orientationFi=1+round(or4*numOrient);                 % 8 is the number of orientations for the filters
orientationFi(orientationFi==(numOrient+1))=1;          % last 'one' is needed to keep correct balance

if size(dataClass,3)==1
    orientationTheta=[];
else
    orientationTheta=zeros(size(orientationFi));
    %[difCols,difRows,difLevs]=gradient(dataClass);
    difLevs=sign(diff(dataClass,1,3));
    difLevs(:,:,end+1)=zeros(rPyr,cPyr,1);   
    rhoAxis=([difCols zeros(rPyr,1,lPyr)].^2+[difRows;zeros(1,cPyr,lPyr)].^2);
    %rhoAxis=([difCols ].^2+[difRows].^2);
    orientationTheta((difLevs>0))=2;
    orientationTheta((difLevs<0))=4;
    orientationTheta(rhoAxis==0)=1;   
    orientationTheta(difLevs==0)=3;
end

% % prodDif=difCols.*difRows;
% % orientationFi(prodDif==0)=3;
% % orientationFi(prodDif>0)=2;
% % orientationFi(prodDif<0)=4;
% % orientationFi(difCols==0)=1; 

%----------------------------------------------------
%------------ subfunction  butterFilter2  ------------
%------------ filter the relevant data witha butterfly filter
%without rearranging, just convolution 
%--------------------------------------------------- 
function    [dataClass] =  butterFilter2 (qtdata,dataClass,boundary,orientationFi,orientationTheta,rPyr,cPyr,lPyr)
%------ input :     qtdata data to be classified
%------             means - prototype of the classes
%------             dataClass las classification
%------             boundary - actual pixels to be filtered 
%------ output:     smooth classification
%------ implicit variables : butterf - the filter to be used for every boundary pixesl    
%------ i                           :  do not pad   the data
%----- first check if there IS a boundary, if not there is only one class
%----- and there is no need of refinement
if sum(sum(sum(boundary)))>0 
    delta=1e-10;
    global dampParam numFeats numClasses butterf levsPyr means algor alpha1 beta1
    %----- dimensions of the variables ------------------------------------------- 
    isData3D=size(dataClass,3);             sizButterf=size(butterf);               
    xSizeButterf=sizButterf(1);             numNeighbours=prod(sizButterf(1:3));
    boundaryCol=boundary(:);                boundaryPixels=sum(boundaryCol~=0);
    %----- if a point is not part of the boundary, take its class mean value (prototype)
    if isData3D>1        reppattern=[1 1 1 numFeats];      reppattern2=[rPyr cPyr lPyr 1];
    else        reppattern=[1 1 numFeats];        reppattern2=[rPyr cPyr  1];
    end
    boundaryMultiD=repmat(boundary,reppattern);
    dataClassMultiD=repmat(dataClass,reppattern);
    qtdata(boundaryMultiD==0)=0;
    for cClass=1:numClasses
        qtdata=qtdata+(boundaryMultiD==0).*(dataClassMultiD==cClass).*repmat(reshape(means(:,cClass),reppattern),reppattern2);
    end
    
    onesNumClasses=ones(1,numClasses);
    maxMeans=(max(means')')*onesNumClasses;
    %--- do not pad, do not check the edges
    numButterf=5;
    for rCoun=3:rPyr-2
        for cCoun=3:cPyr-2
            for lCoun=3:lPyr-2
                if boundary(rCoun,cCoun,lCoun)>0
                    side1Butterf=zeros(125,numFeats);side2Butterf=zeros(125,numFeats);
                    cDims=cCoun-2:cCoun+2;
                    rDims=rCoun-2:rCoun+2;
                    lDims=lCoun-2:lCoun+2;
                    %------ The butterfly filters must be placed correctly
                    filtB=butterf(:,:,:,orientationFi(rCoun,cCoun,lCoun),orientationTheta(rCoun,cCoun,lCoun));  %pick the correct 3D butterfly
                    for cFeats=1:numFeats
                        side1Butterf(:,cFeats)=reshape(qtdata(rDims,cDims,lDims,cFeats).*(filtB==1),5*5*5,1);
                        side2Butterf(:,cFeats)=reshape(qtdata(rDims,cDims,lDims,cFeats).*(filtB==3),5*5*5,1);
                    end
                    currentDataClass1=reshape(dataClass(rDims,cDims,lDims),125,1);
                    currentDataClass2=currentDataClass1;
                    filtB2=reshape(filtB,125,1);
                    currentDataClass1(filtB2~=1)=[];
                    currentDataClass2(filtB2~=3)=[];
                    side1Butterf(filtB2~=1,:)=[];
                    side2Butterf(filtB2~=3,:)=[];
                    
                    [mod1]=mode(currentDataClass1);
                    [mod3]=mode(currentDataClass2);
                    if mod1~=mod3
                        side1Butterf(reshape((currentDataClass1~=mod1),length(currentDataClass1),1),:)=[];
                        side2Butterf(reshape((currentDataClass2~=mod3),length(currentDataClass2),1),:)=[];
                        miu1=mean(side1Butterf);
                        miu2=mean(side2Butterf);
                        miusq1=sum(side1Butterf.^2)/length(side1Butterf);
                        miusq2=sum(side2Butterf.^2)/length(side2Butterf);
                        var1=miusq1-miu1.^2;
                        var3=miusq2-miu2.^2;
                        damp0=((miu1-miu2).^2)./(var1+var3+delta);
                        dampParam=1./(1+exp(alpha1*(beta1-1*damp0)));
                        r1=dampParam'.* squeeze(qtdata(rCoun,cCoun,lCoun,:))+ (1-dampParam').*miu1';
                        r3=dampParam'.* squeeze(qtdata(rCoun,cCoun,lCoun,:)) + (1-dampParam').*miu2';
                        %----compare ONLY between two classes, those which  are most common on each butterfly
                        %for normalised distances divide  by the variance of the POINTS 
                        %%%----- un-normalised distance ---------------
                        dist11=sum(((r1-means(:,mod1)).^2));
                        dist13=sum(((r1-means(:,mod3)).^2));
                        dist31=sum(((r3-means(:,mod1)).^2));
                        dist33=sum(((r3-means(:,mod3)).^2));
                        %           dist11=sum(((r11-means(:,mod1)').^2)./varPoints1,2);
                        %             dist13=sum(((r11-means(:,mod3)').^2)./varPoints1,2);
                        %             dist31=sum(((r33-means(:,mod1)').^2)./varPoints3,2);
                        %             dist33=sum(((r33-means(:,mod3)').^2)./varPoints3,2);
                        [dataMinim1,pixClass1]=min([dist11';dist13';dist31';dist33']); 
                        %------ class assignment --------------------
                        if pixClass1==1 dataClass(rCoun,cCoun,lCoun)=mod1; end
                        if pixClass1==3 dataClass(rCoun,cCoun,lCoun)=mod1; end
                        if pixClass1==2 dataClass(rCoun,cCoun,lCoun)=mod3; end
                        if pixClass1==4 dataClass(rCoun,cCoun,lCoun)=mod3; end
                    else
                        dataClass(rCoun,cCoun,lCoun)=mod3;
                    end
                end
            end
        end
    end
end




%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------
