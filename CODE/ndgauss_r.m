function [ndgauss,sigmax, sigmay,sigmaz]=ndgauss_r(xdim,ydim,zdim,sigmax,sigmay,sigmaz)
% NDGAUSS_R produces an N-dimensional gaussian function (N=1,2,3)
%-----------------------------------------------------------------
% ******         N Dimensional Gauss Function        *******
%-----------------------------------------------------------------
%------  Author :       Constantino Carlos Reyes-Aldasoro
%------                 PHD Student at the University of Warwick
%------  Supervisor :   Abhir Bhalerao
%------  18 October 2001
%-----------------------------------------------------------------
% input   dimensions x,y,z  sigma values x,y,z
% output  n-dimensional gaussian function
%-----------------------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1 help ndgauss_r;  ndgauss=[]; return; end;

%-----------------------------------------------------------------
%------ cases of input:                  -------------------------
%-----------------------------------------------------------------
% 1 only dimensions are specified, then set sigma so that borders 
%        are 1% of central value
%        sigma is set after miu is calculated
% 2 all input arguments are provided



%-----------------------------------------------------------------
%------ dimensions can be input vectors  -------------------------
if nargin==1
   [wRow,wCol,wLev]=size(xdim);
   if wCol==3       %------ 3 D is required, 
      zdim=xdim(3);      ydim=xdim(2);      xdim=xdim(1);
   elseif wCol==2   %------ 2 D is required, set other =1
      ydim=xdim(2);      xdim=xdim(1);      zdim=1;
   elseif wCol==1   %------ 1 D is required, set others =1
      ydim=1;      zdim=1;
   end;
end;

%-----------------------------------------------------------------
%------ only two Dimensions are required ------ set other to one
if nargin==2     zdim=1;   end;

%-----------------------------------------------------------------
%----- x, y, z dimensions of the filter --------------------------
filter.x=1:ceil(xdim);
filter.y=1:ceil(ydim);
filter.z=1:ceil(zdim);
filter.data=zeros(ceil(xdim),ceil(ydim),ceil(zdim));

% mean values for each dimension 
miux=sum(filter.x)/length(filter.x);
miuy=sum(filter.y)/length(filter.y);
miuz=sum(filter.z)/length(filter.z);
%-----------------------------------------------------------------
%set value for 0.1% --> sqrt(2*log(0.001)) = 3.7169
%set value for 1% --> sqrt(2*log(0.01)) = 3.0349
%set value for 10% --> sqrt(2*log(0.1)) = 2.1460
%set value for 50% --> sqrt(2*log(0.5)) = 1.1774
%set value for 70.71% --> sqrt(2*log(0.7071)) = 0.8326   i.e.  50% at the corners


%sigmaDivisor=3.7169;
%sigmaDivisor=3.0349;
%sigmaDivisor=2.1460;
sigmaDivisor=1.1774;
%sigmaDivisor=0.8326;

%------ if no sigmas are provided the set them to a level -----
if nargin < 6
        sigmax=(miux-1)/sigmaDivisor;
        sigmay=(miuy-1)/sigmaDivisor;
        sigmaz=(miuz-1)/sigmaDivisor;
end

%------ sigma must be greater than zero
sigmax=max(sigmax,0.000001);
sigmay=max(sigmay,0.000001);
sigmaz=max(sigmaz,0.000001);

filter.x2=exp(-((filter.x-miux).^2)/2/sigmax/sigmax);
filter.y2=exp(-((filter.y-miuy).^2)/2/sigmay/sigmay);
filter.z2=exp(-((filter.z-miuz).^2)/2/sigmaz/sigmaz);
%filter.xy=exp(-((filter.x-miux).*(filter.y-miuy))/2/sigmax/sigmay);
%filter.yx=exp(-((filter.y-miuy).*(filter.x-miux))/2/sigmax/sigmay);

twodfilter=filter.x2'*filter.y2;
%twodfilter2=filter.xy'*filter.yx;
if ceil(zdim)>1
    for ii=1:ceil(zdim);
        filter.data(:,:,ii)=twodfilter.*filter.z2(ii);
    end;
    ndgauss=filter.data;
else
    ndgauss=twodfilter;
end;

%------ normalising the output ------------------

ndgauss=ndgauss/max(max(max(ndgauss)));