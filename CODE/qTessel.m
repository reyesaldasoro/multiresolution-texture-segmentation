function [q1,q2,q3,q4]=qTessel(data)
%function [q1,q2,q3,q4]=qTessel(data)
%----------------------------------------------------------------------
%****** Quadrant Tessellation Operator  *************------------------
%----------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro------------------
%------             PHD     the University of Warwick------------------
%------  Supervisor :   Abhir Bhalerao    -----------------------------
%------  10 November 2001 ---------------------------------------------
%----------------------------------------------------------------------
%------ Tesselation into Quadrants                  -------------------
%------ Specially to use with MFT algorithm  --------------------------
%------ Based on Wilson and Spann, Image Segmenetation and uncertainty-
%----------------------------------------------------------------------
%------ Input :  Data to be tessellated  must be even numbers ---------
%------ Output : Tessellated data in sections                 ---------  
%------         (1) q1 NW quadrant                            ---------
%------         (2) q2 NE quadrant      q1  |  q2             ---------
%------         (3) q3 SW quadrant      ---------             ---------
%------         (4) q4 SE quadrant      q3  |  q4             ---------
%----------------------------------------------------------------------

%------ no input data is received, error ------------------------------
if nargin<1
   %help qTessel; 
   q1=[];q2=[];q3=[];q4=[];return;
end;

%------ Revise dimensions of data -------------------------------------
[rows,cols,levels,another]=size(data);

%------ !!!!restriction!!! the data sizes/2 must be even.  ------------
if ((rows/2~=floor(rows/2))|(cols/2~=floor(cols/2)))
    warning('Data dimensions must be even');
    centre=[]; surround=[]; help qTessel; return;
end
%------ !!!!restriction!!! the data sizes must at least 8  ------------
%if (rows<8)|(cols<8)
%    warning('Data dimensions too small for operator');
%    centre=[]; surround=[]; help qTessel; return;
%end

%----------------------------------------------------------------------
rowCentr=(rows/2);
colCentr=(cols/2);
q1=data(1:rowCentr,1:colCentr);
q2=data(1:rowCentr,1+colCentr:cols);
q3=data(1+rowCentr:rows,1:colCentr);
q4=data(1+rowCentr:rows,1+colCentr:cols);
