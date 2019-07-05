function [centre,surround,zerCentr]=cTessel(data,dWide)
%function [centre,surround,zerCentr]=cTessel(data,dWide)
%----------------------------------------------------------------------
%****** Centre-surround Tessellation Operator  ******------------------
%----------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro------------------
%------             PHD     the University of Warwick------------------
%------  Supervisor :   Abhir Bhalerao    -----------------------------
%------  10 November 2001 ---------------------------------------------
%----------------------------------------------------------------------
%------ Tesselation into sections Centre - surround -------------------
%------ Specially to use with MFT algorithm  --------------------------
%------ Based on Wilson and Spann, Image Segmenetation and uncertainty-
%----------------------------------------------------------------------
%------ Input :  Data to be tessellated  must be even numbers ---------
%------          dWide determines the wide of centre,1/8 default  -----
%------          but no bigger than 1/2                       ---------
%------ Output : Tessellated data in sections                 ---------
%------         (1) Centre   1/dWide from the central division    -----
%------         (2) Surround all the rest                     ---------
%------         (3) Centre surrounded by zeros                ---------
%----------------------------------------------------------------------

%------ no input data is received, error ------------------------------
if nargin<1   %help cTessel; 
    centre=[]; surround=[];return; 
end;
%------ no input dWide is received,set to 1/8 -------------------------
if nargin==1    dWide=1/8;     end;

%------ Revise dimensions of data -------------------------------------
[rows,cols,levs,another]=size(data);

%------ !!!!restriction!!! the data sizes must be even.  --------------
if ((rows/2~=floor(rows/2))|(cols/2~=floor(cols/2))|(levs/2~=floor(levs/2)&levs~=1))
    warning('Data dimensions must be even');
    centre=[]; surround=[]; help cTessel; return;
end
%------ !!!!restriction!!! the data sizes must be at least 8  ------------
if (rows<4)|(cols<4)|(levs<4&levs~=1)
    warning('Data dimensions too small for operator');
    centre=[]; surround=[]; help cTessel; return;
end
%------ !!!!restriction!!! dWide should be <= than 1/2  ------------
if (dWide>1/2)
    warning('dWide dimensions too big for operator');
    centre=[]; surround=[]; help cTessel; return;
end

%----------------------------------------------------------------------
%----------- Algorithm to remove the Centre  --------------------------
%----------------------------------------------------------------------

%---------Central coordinates ------------------------
rowCentr=(rows/2); colCentr=(cols/2); 

%---------centre width ------------------------
rowCentrWide=floor(rows*dWide);
colCentrWide=floor(cols*dWide);

%---------centre position ------------------------
rowPos=rowCentr+1-rowCentrWide:rowCentr+rowCentrWide;
colPos=colCentr+1-colCentrWide:colCentr+colCentrWide;

%------ special case for levesl because it can be a 2D or 3D data
if levs==1
    levPos=1;
else
    levCentr=(levs/2);
    levCentrWide=max(1,floor(levs*dWide));
    levPos=levCentr+1-levCentrWide:levCentr+levCentrWide;
end

%-------  this is the central data
centre=data(rowPos,colPos,levPos);

data(rowPos,colPos,levPos)=0;
surround=data;
zerCentr=zeros(rows,cols,levs);
zerCentr(rowPos,colPos,levPos)=centre;







