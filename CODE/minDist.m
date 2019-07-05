function [dataClass]=minDist(dataToClassify,prototypes)  
%function [dataClass]=minDist(dataToClassify,prototypes)  
%-----------------------------------------------------------------
%******   minDist : determines minimum distance and assign to a class******-------------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  12 November  2003 ---------------------------------------
%-----------------------------------------------------------------
%------ input  :  DATA / MEANS                       -----  
%------           dataToClassify  is always        2D  [number of points, number of features]
%------           means/prototypes  2D  [number of features,number of classes]
%------ output :  The data classified into the classes predefined   
%-----------------------------------------------------------------
%------  this is a fundamental block for kmeans / kNN / LBG and other classifiers
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<2   help minDist;    dataClass=[];  return; end;
% %------------------------------------------------------------------------
% %------ data has been rearranged, everything into column form  for 2D /3D cases ---
% %------------------------------------------------------------------------
%if ~exist('noNorm') noNorm=0; end
[RCL,numFeats]=size(dataToClassify);
%------ prototypes ------
[numFeats2,numClasses]=size(prototypes);
if (numFeats~=numFeats2)
    warning('dimensions of prototypes are not consistent with dimensions of data');    dataClass=[]; return; 
end
dataClass=zeros(RCL,1);
deltaV=1e-10;
%------ stdData will keep standard deviation of the data per level
%------ it is used to normalise the data and the distances
stdData=deltaV+std(dataToClassify);

    euclDistance=zeros(RCL,numClasses);
    onesDist=ones(RCL,1);   
    %------- these is the normalising stage, to avoid normalising comment these two lines
   % if noNorm==0
        dataToClassify=dataToClassify./(onesDist*stdData);                  %-----NORMALISE WITH STD
        prototypes=prototypes./(stdData'*ones(1,numClasses));     %-----NORMALISE WITH STD
   % end
    %-------  euclidean distancecompute relative distance 
    for counter2=1:numClasses
%         if size(dataToClassify,1)>4e6
%             for counterFeats=1:numFeats
%                 distanceToPrototype(:,counterFeats)=dataToClassify(:,counterFeats)-onesDist*(prototypes(counterFeats,counter2)');
%                 distance2(:,counterFeats)=distanceToPrototype(:,counterFeats).^2;         %---- For Euclidean
%                 %euclDistance(:,counter2)=sum(distance2,2);
%             end
%         else
            distanceToPrototype=dataToClassify-onesDist*(prototypes(:,counter2)');
            distance2=distanceToPrototype.^2;         %---- For Euclidean
            %distance2=abs(distanceToPrototype);       %---- For City Block
            %euclDistance(:,counter2)=sum(distance2,2);
%        end
        euclDistance(:,counter2)=sum(distance2,2);
    end

    %------ find the minimum and thus to which class it belongs
    %------ and on the same fly calculate the error up to the feature set considered        
    [dataMinim,dataClass]=(min(euclDistance'));
    dataClass=dataClass';
