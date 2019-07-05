function [data]=qDeTessel(q1,q2,q3,q4)
%function [data]=qDeTessel(q1,q2,q3,q4)
%----------------------------------------------------------------------
%****** Quadrant DE - Tessellation Operator  *************------------------
%----------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro------------------
%------             PHD     the University of Warwick------------------
%------  Supervisor :   Abhir Bhalerao    -----------------------------
%------  16 November 2001 ---------------------------------------------
%----------------------------------------------------------------------
%------ De Tesselation into Quadrants                  -------------------
%------ Specially to use with MFT algorithm  --------------------------
%------ Based on Wilson and Spann, Image Segmenetation and uncertainty-
%----------------------------------------------------------------------
%------ Input : Tessellated data in sections                 ---------  
%------         (1) q1 NW quadrant                            ---------
%------         (2) q2 NE quadrant      q1  |  q2             ---------
%------         (3) q3 SW quadrant      ---------             ---------
%------         (4) q4 SE quadrant      q3  |  q4             ---------
%------ Output : Data reassembled from tessellation           ---------
%----------------------------------------------------------------------

%------ no input data is received, error ------------------------------
if nargin<1
   %help qDeTessel; 
   data=[];return;
end;

%------ Revise dimensions of data -------------------------------------
switch nargin
case 4
    [q4rows,q4cols,q4levels]=size(q4);
    [q3rows,q3cols,q3levels]=size(q3);
    [q2rows,q2cols,q2levels]=size(q2);
    [q1rows,q1cols,q1levels]=size(q1);
    if q1==0
        q1=zeros(max([q4rows q3rows q2rows q1rows]),max([q4cols q3cols q2cols q1cols]));
    end
    if q2==0
        q2=zeros(max([q4rows q3rows q2rows q1rows]),max([q4cols q3cols q2cols q1cols]));
    end
    if q3==0
        q3=zeros(max([q4rows q3rows q2rows q1rows]),max([q4cols q3cols q2cols q1cols]));
    end
    if q4==0
        q4=zeros(max([q4rows q3rows q2rows q1rows]),max([q4cols q3cols q2cols q1cols]));
    end
   
case 3
    q4=zeros(size(q1));
case 2
    q3=zeros(size(q1));
    q4=zeros(size(q1));
case 1
    q2=zeros(size(q1));
    q3=zeros(size(q1));
    q4=zeros(size(q1));
otherwise
    %help qDeTessel; 
    data=[];return;
end;


%----------------------------------------------------------------------
try
    data=[q1 q2;q3 q4];
catch
    warning ('Data dimensions do not agree');
    %help qDeTessel; 
    data=[];return;
end
