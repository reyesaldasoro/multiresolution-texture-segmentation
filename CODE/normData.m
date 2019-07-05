function [dataNorm]=normData(dataToNormalise)  
%function  [dataNorm]=normData(dataToNormalise)  
%-----------------------------------------------------------------
%******  normalises data by dividing over the std 
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  27 November  2003 ---------------------------------------
%-----------------------------------------------------------------
%------ input  :  DATA                      -----  
%------           dataToNormalise  is always        2D  [number of points x number of features]
%------ output :  The data normalised to have a new STD of 1
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin~=1   help normData;    dataNorm=[];  return; end;
% %------------------------------------------------------------------------
% %------ data has been rearranged, everything into column form  for 2D /3D cases ---
% %------------------------------------------------------------------------

[RCL,numFeats]=size(dataToNormalise);
deltaV=1e-10;
%------ stdData will keep standard deviation of the data per level
%------ it is used to normalise the data and the distances
stdData=deltaV+std(dataToNormalise);
onesDist=ones(RCL,1);   
    %------- these is the normalising stage, to avoid normalising comment these two lines
    dataNorm=dataToNormalise./(onesDist*stdData);                  %-----NORMALISE WITH STD
