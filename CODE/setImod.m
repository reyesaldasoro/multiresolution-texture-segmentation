function Imodn=setImod(Imod,n)
% function Imodn=setImod(Imod,n)
%----- function that calculates the final Imod neccessary
%----- a recursive function is used
%----------------------------------------------------------------
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------   November  2001 ----------------------------
%----------------------------------------------------
%------ input : Imod matrix with dimensions defined -
%------         for the uniform or gaussian cases  --
%------         n  necessarily even number for uniform

%------ no input data is received, error -------------------------
if nargin<1
   help setImod;
   Imodn=[];
   return;
end;

if (n/2)~=floor(n/2)
   help setImod;
   Imodn=[];
   return;
end;


%----- for the uniform reduction/expansion a Imod would
%----- be a Identity matrix duplicated the rows or cols
%----- the relationship between lines and cols is *2*
%----- the Gaussian Imod would replicate the kernel function
%----- and therefore relationship will be different 1-5, 2-7, 3-9...

[a,b]=size(Imod);
I1=max(a,b);
I2=min(a,b);

if (I1/I2)==2         %this would imply the uniform case
  if size(Imod,1)==n
      Imodn=Imod;
  else
      Imodn=setImod([Imod zeros(size(Imod));zeros(size(Imod)) Imod],n);
  end;
elseif (I1==(3+2*I2))  %this would imply the gaussian case REDUCTION
  if I2==n
      Imodn=Imod;
  else
      if (a>b)
         Imodn=setImod([[Imod;zeros(2*b,b)] [zeros(2*b,b);Imod]],n);
      else
         Imodn=setImod([[Imod zeros(a,2*a)]; [zeros(a,2*a) Imod]],n);
      end;
  end;
else
  if b/2==n
      Imodn=Imod;
  else
      Imodn=setImod([[Imod;zeros(b/2,b)] [zeros(b/2,b);Imod]],n);
  end;

end;
