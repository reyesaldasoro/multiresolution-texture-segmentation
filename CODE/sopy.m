function [d]=sopy(data,levelsDown,inversefft)
%function [d]=sopy(data,levelsDown,inversefft)
%-----------------------------------------------------------------
%******  Second Orientation PYramid            ******-------------
%------  sopy transforms data into the Fourier Domain and then
%------  filters it with a Second Orientation Pyramid tessellation
%------  with truncated Gaussians in different frequency-orientation positions
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  5 March 2002    -----------------------------------------
%-----------------------------------------------------------------
%------ input  :  The data  upon which a second orientation  -----
%------           Pyramid will be obtained                   -----
%------           data will be received as spatial domain    -----
%------           optional param: num levels down            -----
%------ output :  The different levels of the Pyramid        -----
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1   help sopy;    d=[];   return; end;
[rows,cols,levs]=size(data);

%if nargin==1 levelsDown=floor(log2(rows))-2; end;
if ~exist('levelsDown') levelsDown=1; end
if ~exist('inversefft') inversefft=1; end

%-------------frequency pyramid analisis zooming
if levs>1 warning('caution, data is 3d'); d=sopy3d(data,levelsDown); return; end   %------ call sopy3d in case data is 3d

%dataf=fftshift(fft2(fftshift(data)));  
dataf=fftshift(fft2((data)));  

%----- This shift of the data on the fourier domain is important in the
%-----    cases when the energy is concentrated on the line n/2+1
dataf=dataf(rows:-1:1,cols:-1:1,levs:-1:1);



%-------------first level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c1,s2_7]=tessel6(dataf);               %---- 512   center 256  BP 128
d(:,:,1:7)=deTessel6(c1,s2_7,rows,inversefft);
levelsDown=levelsDown-1;
%-------------second level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c8,s9_14]=tessel6(c1);                 %---- 256   center 128  BP 64
d(:,:,8:14)=deTessel6(c8,s9_14,rows,inversefft);
levelsDown=levelsDown-1;
%-------------third level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c15,s16_21]=tessel6(c8);               %---- 128   center 64   BP 32
d(:,:,15:21)=deTessel6(c15,s16_21,rows,inversefft);
levelsDown=levelsDown-1;
%-------------fourth level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c22,s23_28]=tessel6(c15);              %----  64   center 32   BP 16
d(:,:,22:28)=deTessel6(c22,s23_28,rows,inversefft);
levelsDown=levelsDown-1;
%-------------fifth level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c29,s30_35]=tessel6(c22);              %----  32   center 16   BP 8
d(:,:,29:35)=deTessel6(c29,s30_35,rows,inversefft);
levelsDown=levelsDown-1;
%-------------sixth level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c36,s37_42]=tessel6(c29);              %----  16   center 8    BP 4
d(:,:,36:42)=deTessel6(c36,s37_42,rows,inversefft);
levelsDown=levelsDown-1;
%-------------seventh level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c43,s44_49]=tessel6(c36);              %----   8   center 4    BP 2
d(:,:,43:49)=deTessel6(c43,s44_49,rows,inversefft);
levelsDown=levelsDown-1;
%-------------eightth level ------------------------
if levelsDown==0 d=d(rows:-1:1,cols:-1:1,:); return; end;
[c50,s51_56]=tessel6(c43);              %----   4   center 2    BP 1   
d(:,:,50:56)=deTessel6(c50,s51_56,rows,inversefft);

%------------------------------------------------------
%---------  Subfunctions         ----------------------
%------------------------------------------------------

function [q1,q2]=tessel6(dataf)
%----------------------------------------------------------------------
%****** 7 section (C+6 surroundings) Tessellation Operator  ******-----
%----------------------------------------------------------------------

[q1,s]=cTessel(dataf,1/4);     %-----------Tessellation process
[s1,s2,s3,s4]=qTessel(s);
[qt3,qt4,qt2,qx]=qTessel(s1);
[qt5,qt6,qx,qt7]=qTessel(s2);

q2(:,:,1)=qt2;  q2(:,:,2)=qt3;  q2(:,:,3)=qt4;  
q2(:,:,4)=qt5;  q2(:,:,5)=qt6;  q2(:,:,6)=qt7;


function [d]=deTessel6(q1,q2,n,inversefft)
%----------------------------------------------------------------------
%----- n determines the final size of each segmentation
d1 =cDeTessel(ndgauss_r(size(q1)).*q1,1/4);
gFilter=ndgauss_r(size(q2(:,:,1)));
d2 =qDeTessel(qDeTessel(0,0,gFilter.*q2(:,:,1),0),0,0,0);
d3 =qDeTessel(qDeTessel(gFilter.*q2(:,:,2),0,0,0),0,0,0);
d4 =qDeTessel(qDeTessel(0,gFilter.*q2(:,:,3),0,0),0,0,0);
d5 =qDeTessel(0,qDeTessel(gFilter.*q2(:,:,4),0,0,0),0,0);
d6 =qDeTessel(0,qDeTessel(0,gFilter.*q2(:,:,5),0,0),0,0);
d7 =qDeTessel(0,qDeTessel(0,0,0,gFilter.*q2(:,:,6)),0,0);

%----- if the segmentations are not yet the size, then surround with zeros
[a,b]=size(d1);
if a<n
   d1=cDeTessel(d1,a/n/2);   d2=cDeTessel(d2,a/n/2);   d3=cDeTessel(d3,a/n/2);
   d4=cDeTessel(d4,a/n/2);   d5=cDeTessel(d5,a/n/2);   d6=cDeTessel(d6,a/n/2);
   d7=cDeTessel(d7,a/n/2);
end

%------ finally inverse transform if needed!!!!
if inversefft==1
   d(:,:,1)=abs((ifftn(d1)));   d(:,:,2)=abs((ifftn(d2)));   d(:,:,3)=abs((ifftn(d3)));
   d(:,:,4)=abs((ifftn(d4)));   d(:,:,5)=abs((ifftn(d5)));   d(:,:,6)=abs((ifftn(d6)));
   d(:,:,7)=abs((ifftn(d7)));
%    d(:,:,1)=abs(fftshift(ifftn(d1)));   d(:,:,2)=abs(fftshift(ifftn(d2)));   d(:,:,3)=abs(fftshift(ifftn(d3)));
%    d(:,:,4)=abs(fftshift(ifftn(d4)));   d(:,:,5)=abs(fftshift(ifftn(d5)));   d(:,:,6)=abs(fftshift(ifftn(d6)));
%    d(:,:,7)=abs(fftshift(ifftn(d7)));
else
   d(:,:,1)=(ifftn(d1));   d(:,:,2)=(ifftn(d2));   d(:,:,3)=(ifftn(d3));
   d(:,:,4)=(ifftn(d4));   d(:,:,5)=(ifftn(d5));   d(:,:,6)=(ifftn(d6));
   d(:,:,7)=(ifftn(d7));
%    d(:,:,1)=fftshift(ifftn(d1));   d(:,:,2)=fftshift(ifftn(d2));   d(:,:,3)=fftshift(ifftn(d3));
%    d(:,:,4)=fftshift(ifftn(d4));   d(:,:,5)=fftshift(ifftn(d5));   d(:,:,6)=fftshift(ifftn(d6));
%    d(:,:,7)=fftshift(ifftn(d7));
end

function [e]=normalise(dd,data)
%----------------------------------------------------------------------
%----- normalises the d's w.r.t. data and substracts each from data

totData=sum(sum(data));
e(:,:,1)=(data-abs(dd(:,:,1))*totData/(sum(sum(abs(dd(:,:,1))))));
e(:,:,2)=(data-abs(dd(:,:,2))*totData/(sum(sum(abs(dd(:,:,2))))));
e(:,:,3)=(data-abs(dd(:,:,3))*totData/(sum(sum(abs(dd(:,:,3))))));
e(:,:,4)=(data-abs(dd(:,:,4))*totData/(sum(sum(abs(dd(:,:,4))))));
e(:,:,5)=(data-abs(dd(:,:,5))*totData/(sum(sum(abs(dd(:,:,5))))));
e(:,:,6)=(data-abs(dd(:,:,6))*totData/(sum(sum(abs(dd(:,:,6))))));
e(:,:,7)=(data-abs(dd(:,:,7))*totData/(sum(sum(abs(dd(:,:,7))))));






