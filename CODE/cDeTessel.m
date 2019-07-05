function [data]=cDeTessel(c,dWide)
%function [data]=cDeTessel(c,dWide)
%----------------------------------------------------------------------
%****** Centre DE - Tessellation Operator  *************------------------
%----------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro------------------
%------             PHD     the University of Warwick------------------
%------  Supervisor :   Abhir Bhalerao    -----------------------------
%------  16 November 2001 ---------------------------------------------
%----------------------------------------------------------------------
%------ De Tesselation in centre/surround           -------------------
%------ Specially to use with MFT algorithm  --------------------------
%------ Based on Wilson and Spann, Image Segmenetation and uncertainty-
%----------------------------------------------------------------------
%------ Input : Tessellated data in sections                 ---------  
%------         (1) Centre   1/dWide from the central division    -----
%------         (2) Surround all the rest                     ---------
%------ Output : Data reassembled from tessellation           ---------
%----------------------------------------------------------------------

%------ no input data is received, error ------------------------------
if nargin<1  %help cDeTessel; 
    data=[];return; end;
%------ no input dWide is received,set to 1/8 -------------------------
if nargin==1    dWide=1/4;     end;

%------ Revise dimensions of data -------------------------------------
[rows,cols,levs,another]=size(c);
%------ determine size of the original image from which the centre was extracted
imRows=rows/dWide/2;
imCols=cols/dWide/2;
if levs==1
    imLevs=1;
    data=zeros(imRows,imCols,imLevs);
    %------ Place data in the centre
    data(imRows/2+1-rows/2:imRows/2+rows/2,imCols/2+1-cols/2:imCols/2+cols/2)=c;
else
    imLevs=levs/dWide/2;
    data=zeros(imRows,imCols,imLevs);
    %------ Place data in the centre
    data(imRows/2+1-rows/2:imRows/2+rows/2,imCols/2+1-cols/2:imCols/2+cols/2,imLevs/2+1-levs/2:imLevs/2+levs/2)=c;
end

%size(data);

%----------------------------------------------------------------------
