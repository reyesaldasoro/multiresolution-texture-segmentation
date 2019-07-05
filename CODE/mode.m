function [modeVal,modeVal2]=mode(datain,modeDim)
%------  function [modeVal,modeVal2]=mode(datain,modeDim)
%------  Function that gets the mode of a certain vector or matrix
%------  if modeDims is not present then the mode is as a whole, otherwise
%-----   take it by columns
%-----------------------------------------------------------------
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  6 November 2001 ----------------------------
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1
  help mode; modeVal=[]; return;
end;

[rows,cols,levels]=size(datain);
mind=min(min(min(datain)));
maxd=max(max(max(datain)));

if sum(sum(datain~=floor(datain)))~=0    stepAxis=0.1;
else    stepAxis=1;
end
%------  the axis vector which is required by histc and will indicate the mode
xAxisVec=mind:stepAxis:maxd;

if ~exist('modeDim')
    datain=reshape(datain,1,cols*rows*levels);
end

%------  now that is a row vector, use HISTC to count the values and find maximum
[elemCount]=histc(datain,xAxisVec);
[modeCount,modeIndex]=max(elemCount);
modeVal=xAxisVec(modeIndex);

%------- eliminate mode values and find second mode
if ~exist('modeDim')
    elemCount(modeIndex)=0;
else
    for cCol=1:cols elemCount(modeIndex(cCol),cCol)=0; end
end
[modeCount,modeIndex]=max(elemCount);
modeVal2=xAxisVec(modeIndex);


