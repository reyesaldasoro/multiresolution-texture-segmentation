function [centres]=LBG (data,numCentres)
%function [centres]=LBG (data,numCentres)
%-----------------------------------------------------------------
%******  Linde Buzo Gray vector quantising algorothm ******-------------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  1 March 2004    -----------------------------------------
%-----------------------------------------------------------------
%------ input  :  The data for which a number of centres 
%------           will be obtained                           -----
%------           the number of centres                      -----
%------ output :  The different centres                      -----
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

%----- in case data is a series of images 
%----- transform into column data
% [rows,cols,levs,numFeats]=size(data);
% if levs>1
%     if numFeats>1
%         for k=1:numFeats
%             data2(k,:)=reshape(data(:,:,:,k),1,rows*cols*levs);
%         end
%     else
%         for k=1:levs
%             data2(k,:)=reshape(data(:,:,k),1,rows*cols);
%         end
%     end
%     clear data; 
%     data=data2;
%     clear data2;
% end
%data=data';
[data,rows,cols,levs,numFeats,RCL]=arrangeData(data);
epsilon=0.051;
chainL={'ro','bo','ko','go','co','mo'};
%----- initial values -----
[RCL,numFeats]=size(data);
centres(1,:)=(mean(data));
data_Centred=data-repmat(centres,RCL,1);
distortion_old=sum(sum(((data_Centred).^2)));
diffDistortions=distortion_old;
%figure;
%xyclass(data(:,1),data(:,2),'','','k.',2);hold
%xyclass(centres(:,1),centres(:,2),'','','m*',5);

%------ outer cycle of splitting centres
for counterSplits=1:ceil(log2(numCentres))
    %splitting
    lenCentres=size(centres,1);
    for counterC=1:lenCentres
        centres(counterC+lenCentres,:)=   (1-epsilon)*centres(counterC,:);
        centres(counterC,:)=              (1+epsilon)*centres(counterC,:);
    end
    %figure;
    %xyclass(data(:,1),data(:,2),'','','k.',2);hold
    %xyclass(centres(:,1),centres(:,2),'','','m*',5);
    lenCentres=lenCentres*2;
    %second cycle where the new centres move to minimise the distortion	
    while (diffDistortions>epsilon)
        % i assign a class
        for counterC=1:lenCentres
            data_Centred=data-repmat(centres(counterC,:),RCL,1);
            distance(:,counterC)= sum(data_Centred.^2,2);
        end
        
        [minDist,pixClass]=min(distance,[],2);
        % ii Update the codevector (centres)
        for counterC=1:lenCentres
            if sum(pixClass==counterC)>0  %---- if there are no pixels in the class, keep the same position
                data_ClassX=data(pixClass==counterC,:);
                %xyclass(data(pixClass==counterC,1),data(pixClass==counterC,2),'','',chainL{counterC},5)
                centres(counterC,:)=mean(data_ClassX);
            end
        end
        %xyclass(centres(:,1),centres(:,2),'','','g*',5);
        % iv  get new distortions
        distortion_new=0;
        for counterC=1:lenCentres
            data_ClassX=data(pixClass==counterC,:);
            pixInClass=sum(pixClass==counterC);
            data_Centred=data_ClassX-repmat(centres(counterC,:),pixInClass,1);
            distortion_new=distortion_new+sum(sum(data_Centred.^2));
        end
        diffDistortions=(distortion_old-distortion_new)/distortion_old;
        distortion_old=distortion_new;
        %[counterSplits diffDistortions]
    end
    distortion_old=distortion_new;
    diffDistortions=distortion_old;
end
centres=centres';


centres=centres(:,1:numCentres);

