function [dataArranged,rows,cols,levs,numFeats,RCL]=arrangeData(varargin)  
%function [dataArranged]=arrangeData(dataIn)  
%-----------------------------------------------------------------
%****** Rearrangement of data from several formats to 2D [RCL x numFeats]*****-------------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  17 November  2003 ---------------------------------------
%-----------------------------------------------------------------
%------ input  :  DATA                       -----  
%------           Data and Mask can be 2d/3d mask is always required
%------ output :  The data rearranged
%-----------------------------------------------------------------
%------ The parameters of the algorithm are:                  ----
%------ dataSets        : the number of sets as input (read from data)---
%-----------------------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1   help arrangeData;    dataArranged=[];  err=[]; return; end;
% %------------------------------------------------------------------------
% %------ rearrange data, everything into column form  for 2D /3D cases ---
% %------------------------------------------------------------------------
%------ number of data sets received = numParameters
%------  Data input options --------------------------------------
%------     1 - one single (N+1)D matrix of (N)D feature matrices,      %length(varagin)==1
%------     2 - a series of ND feature matrices,                                    %length(varagin)>1
%------ There should be a consistency of the number of sets in the means and the data    
%------ valid options   numFeats=1 , dataSets=1  just one set to classify
%------                 numFeats=n , dataSets=n  several data sets 
%------                 numFeats=n , dataSets=1  several data sets embeded
%------ invalid option  numFeats=1 , dataSets=n  no means for classification
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

if  (length(varargin)==1)&(iscell(varargin{1}))
    dataIn=varargin{1};
    dataSets=length(dataIn);
else
    dataSets=length(varargin);                      
    dataIn=varargin;
end

if dataSets==1                                  %can be just one set or can be multidimensional
    dataIn=varargin{1};
    if iscell(dataIn) dataIn=dataIn{1}; end
    [rows,cols,levs,other]=size(dataIn);
    if other>1                                                  %----- 4D data, a series of 3D cubes
        numFeats=other;
        RCL=rows*cols*levs;
        dataArranged=reshape(dataIn,RCL,numFeats);
    elseif levs>1                                               %----- 3D data a series of 2D slices
        numFeats=levs;
        levs=1;
        RCL=rows*cols;
        dataArranged=reshape(dataIn,RCL,numFeats);
    else                                                        %----- 2D data ... either one image or arranged data
        if rows>(10*cols)                                       %----- if data is roughly square it is an image
            dataArranged=dataIn;                                    % these is in case data receive has already been arranged!
            RCL=rows;numFeats=cols;
            sqrtRCL=sqrt(RCL);
            if floor(sqrtRCL)==sqrtRCL
                rows=sqrtRCL; cols=sqrtRCL; levs=1;
            else
                sqrtRCL=(RCL)^(1/3);
                if floor(sqrtRCL)==sqrtRCL rows=sqrtRCL; cols=sqrtRCL; levs=sqrtRCL;
                else rows=[];cols=[]; levs=[];
                end
            end
        else
            dataArranged=dataIn(:);
            RCL=rows*cols;levs=1;numFeats=1;

        end

    end
    
else                       %several feature data sets
    numFeats=length(dataIn);
    [rows,cols,levs]=size(dataIn{1});
    RCL=rows*cols*levs;
    for counter = 1:numFeats    %----- all the rest of the arguments are considered data unless ones is nxnxM
        dataArranged(:,counter)=reshape(dataIn{counter},RCL,1);
    end
end
