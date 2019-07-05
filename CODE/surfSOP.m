function surfSOP(qt2,qt3,qt4,qt5,qt6,qt7,qt8,qt9,qt10,qt11,qt12,qt13,qt14,numfig)
%function surfSOP(qt1,qt2,qt3,qt4,qt5,qt6,qt7,qt8,qt9,qt10,qt11,qt12,qt13,numfig)
%----------------------------------------------------------------------
%****** Display 2 levels of the SOP results                ******------------------
%----------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro------------------
%------             PHD     the University of Warwick------------------
%------  Supervisor :   Abhir Bhalerao    -----------------------------
%------  12 November 2001 ---------------------------------------------
%----------------------------------------------------------------------

%------ no input data is received, error ------------------------------
if nargin<1   help surfSOP; return;  end;

if nargin==6    qt8=qt1;qt9=qt2;qt10=qt3;qt11=qt4;qt12=qt5;qt13=qt6; end;

if (nargin==1|nargin==13)|nargin==6    numfig=3;  end

if nargin <=2
   %------ data produced directly from MFT.M function --------------------
   %------ Revise dimensions of data -------------------------------------
   dataSOP=qt2;
   if nargin==2        numfig=qt3;    end;
   if iscell(dataSOP)
       qt2=dataSOP{1};       qt3=dataSOP{2};       qt4=dataSOP{3};
       qt5=dataSOP{4};       qt6=dataSOP{5};       qt7=dataSOP{6};
       qt8=dataSOP{7};       qt9=dataSOP{8};       qt10=dataSOP{9};
       qt11=dataSOP{10};       qt12=dataSOP{11};       qt13=dataSOP{12};
       qt14=dataSOP{13};
       levels=13;
       [rows,cols]=size(qt2);
   else
       [rows,cols,levels,another]=size(dataSOP);
       qt2=dataSOP(:,:,1);       qt3=dataSOP(:,:,2);       qt4=dataSOP(:,:,3);
       qt5=dataSOP(:,:,4);       qt6=dataSOP(:,:,5);       qt7=dataSOP(:,:,6);
       qt8=dataSOP(:,:,7);       qt9=dataSOP(:,:,8);       qt10=dataSOP(:,:,9);
       qt11=dataSOP(:,:,10);       qt12=dataSOP(:,:,11);       qt13=dataSOP(:,:,12);
       qt14=dataSOP(:,:,13);
   end     
   
   %------ !!!!restriction!!! the data sizes must be even.  --------------
   if levels~=13        warning('Data dimensions not adequate'); help surfSOP; return;    end;
end;

%------ funcMFT defines the type of plot: absolute, angle, angleunwrapped
%------ only in the cases of complex numbers, otherwise is always abs
if(isreal(qt2))
%    funcSOP='((abs(';
    funcSOP='(((';
else
   funcSOP='((abs(';
%   funcSOP='((real(';
%   funcSOP='((imag(';
%   funcSOP='((angle(';
%    funcSOP='(unwrap2d(angle(';
%    funcSOP='fftshift(fft2(angle(';
%    funcSOP='abs(diff(angle(';
%    funcSOP='abs(diff(angle(';
end
%func2=''')))';
func2=')))';
    
    qt14=eval(strcat(funcSOP,'qt14',func2));
    qt2=eval(strcat(funcSOP,'qt2',func2));
    qt3=eval(strcat(funcSOP,'qt3',func2));
    qt4=eval(strcat(funcSOP,'qt4',func2));
    qt5=eval(strcat(funcSOP,'qt5',func2));
    qt6=eval(strcat(funcSOP,'qt6',func2));
    qt7=eval(strcat(funcSOP,'qt7',func2));
    qt8=eval(strcat(funcSOP,'qt8',func2));
    qt9=eval(strcat(funcSOP,'qt9',func2));
    qt10=eval(strcat(funcSOP,'qt10',func2));
    qt11=eval(strcat(funcSOP,'qt11',func2));
    qt12=eval(strcat(funcSOP,'qt12',func2));
    qt13=eval(strcat(funcSOP,'qt13',func2));

%------  in some cases it is convenient to reduce the size of the plot

stepPlot=2;
% 
% if size(qt14,1)>128 qt14=(((qt14(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt2,1)>128 qt2=(((qt2(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt3,1)>128 qt3=(((qt3(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt4,1)>128 qt4=(((qt4(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt5,1)>128 qt5=(((qt5(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt6,1)>128 qt6=(((qt6(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt7,1)>128 qt7=(((qt7(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt8,1)>128 qt8=(((qt8(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt9,1)>128 qt9=(((qt9(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt10,1)>128 qt10=(((qt10(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt11,1)>128 qt11=(((qt11(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt12,1)>128 qt12=(((qt12(1:stepPlot:end,1:stepPlot:end)))); end;
% if size(qt13,1)>128 qt13=(((qt13(1:stepPlot:end,1:stepPlot:end)))); end;

if size(qt14,1)>128 qt14=reduceu(reduceu(reduceu(reduceu((qt14))))); end;
%if size(qt1,1)>128 qt1=reduceu(reduceu((qt1))); end;
if size(qt2,1)>128 qt2=reduceu(reduceu(reduceu(reduceu((qt2))))); end;
if size(qt3,1)>128 qt3=reduceu(reduceu(reduceu(reduceu((qt3))))); end;
if size(qt4,1)>128 qt4=reduceu(reduceu(reduceu(reduceu((qt4))))); end;
if size(qt5,1)>128 qt5=reduceu(reduceu(reduceu(reduceu((qt5))))); end;
if size(qt6,1)>128 qt6=reduceu(reduceu(reduceu(reduceu((qt6))))); end;
if size(qt7,1)>128 qt7=reduceu(reduceu(reduceu(reduceu((qt7))))); end;
if size(qt8,1)>128 qt8=reduceu(reduceu(reduceu(reduceu((qt8))))); end;
if size(qt9,1)>128 qt9=reduceu(reduceu(reduceu(reduceu((qt9))))); end;
if size(qt10,1)>128 qt10=reduceu(reduceu(reduceu(reduceu((qt10))))); end;
if size(qt11,1)>128 qt11=reduceu(reduceu(reduceu(reduceu((qt11))))); end;
if size(qt12,1)>128 qt12=reduceu(reduceu(reduceu(reduceu((qt12))))); end;
if size(qt13,1)>128 qt13=reduceu(reduceu(reduceu(reduceu((qt13))))); end;


figure(numfig)
%------ external figures first -----------------------
subplot('position',[.01 .2 .23 .38])
imagesc((qt2));  axis off;
subplot('position',[.01 .6 .23 .38])
imagesc((qt3));  axis off;
subplot('position',[.26 .6 .23 .38])
imagesc((qt4));  axis off;
subplot('position',[.51 .6 .23 .38])
imagesc((qt5));  axis off;
subplot('position',[.76 .6 .23 .38])
imagesc((qt6));  axis off;
subplot('position',[.76 .2 .23 .38])
imagesc((qt7));  axis off;

%------ middle figures ------------------------------
if nargin~= 6
    subplot('position',[.26 .2 .11 .18])
    imagesc((qt9));      axis off;
    subplot('position',[.26 .4 .11 .18])
    imagesc((qt10));    axis off;
    subplot('position',[.38 .4 .11 .18])
    imagesc((qt11));    axis off;
    subplot('position',[.51 .4 .11 .18])
    imagesc((qt12));    axis off;
    subplot('position',[.63 .4 .11 .18])
    imagesc((qt13));    axis off;
    subplot('position',[.63 .2 .11 .18])
    imagesc((qt14));    axis off;
    %------ central position ------------------------------
    subplot('position',[.38 .01 .24 .37])
    imagesc((qt8));    axis off;
end;










% %------ external figures first -----------------------
% subplot('position',[.01 .2 .23 .38])
% surfdat((qt8))
% axis off;
% subplot('position',[.01 .6 .23 .38])
% surfdat((qt9));
% axis off;
% subplot('position',[.26 .6 .23 .38])
% surfdat((qt10))
% axis off;
% subplot('position',[.51 .6 .23 .38])
% surfdat((qt11))
% axis off;
% subplot('position',[.76 .6 .23 .38])
% surfdat((qt12))
% axis off;
% subplot('position',[.76 .2 .23 .38])
% surfdat((qt13))
% axis off;
% 
% %------ middle figures ------------------------------
% if nargin~= 6
%     subplot('position',[.26 .2 .11 .18])
%     surfdat((qt2))
%     axis off;
%     subplot('position',[.26 .4 .11 .18])
%     surfdat((qt3))
%     axis off;
%     subplot('position',[.38 .4 .11 .18])
%     surfdat((qt4))
%     axis off;
%     subplot('position',[.51 .4 .11 .18])
%     surfdat((qt5))
%     axis off;
%     subplot('position',[.63 .4 .11 .18])
%     surfdat((qt6))
%     axis off;
%     subplot('position',[.63 .2 .11 .18])
%     surfdat((qt7))
%     axis off;
% 
%     %------ central position ------------------------------
%     subplot('position',[.38 .01 .24 .37])
%     surfdat((qt1))
%     axis off;
