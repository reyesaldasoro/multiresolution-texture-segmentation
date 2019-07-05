function nProb=normProb(x,miux,sigmax)
%function nProb=normProb(x,miux,sigmax)
%-----------------------------------------------------------------
%******   Normal Probability calculation                ***-------
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  3 March 2002     ----------------------------------------
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

[rowm,colm,levm]=size(miux);
[a,b]=size(sigmax);
%----- miux always a column vector
if (colm>1)&(rowm>1)  warning('miu dimensions  not consistent ');  nProb=[];  return; end
if (colm>rowm)&(rowm==1) miux=miux'; rowm=colm; colm=1; end
%----- miux and sigma must have consistent dimensions
if (rowm~=a)  warning('S and miu dimensions  not consistent ');  nProb=[];  return; end
if (b~=a)  warning('S dimensions  not consistent ');  nProb=[];  return; end


if all(size(miux)==1)                  %----- 1D one point or a series of points
    denomin=(sigmax)*((2*pi)^(1/2));
    numerat=-0.5*((x-miux).^2)/sigmax/sigmax;
    nProb=(exp(numerat))/(denomin);
else                                    %----- 2D or more
    %------ for more than 1D the    means have to be Nx1  
    %------                         sigma has to be NxN
    %------                         x     has to be jxkxN
    %------    first rearrange x into a jxk x N matrix  
    [rowX,colX,levX]=size(x);
    if levX==1                          % only 1D or 2D matrices
        if (any([rowX==1 colX==1]))     % either  [p p] or [p;p]
            if colX>rowX x=x'; end
            difmiu=x-miux;
        elseif (any([rowX==2 colX==2])) % either  [2 x colX] or [rowX x 2]
            if colX<rowX x=x'; end
            difmiu(1,:)=x(1,:)-miux(1,1);
            difmiu(2,:)=x(2,:)-miux(2,1);
        elseif (any([rowX>2 colX>2])) % either  [n x colX] or [rowX x n]
            if colX<rowX x=x'; [rowX,colX]=size(x); end
            for counter=1:rowX
                difmiu(counter,:)=x(counter,:)-miux(counter,1);
            end
        end
    else                                % levX~=1  3D Matrices
        if  (rowX==colX)               % [rowX x colX] [rowX x colX] 
            for counter=1:levX
                difmiu(counter,:)=reshape(x(:,:,counter),1,rowX*colX)-miux(counter,1);
            end
            %difmiu(2,:)=reshape(x(:,:,2),1,rowX*colX)-miux(2,1);
            invS=inv(sigmax); detS=det(sigmax);
            denomin=sqrt(detS)*(2*pi)^(a/2);
            numerat1=-0.5*(difmiu)'*invS;
            numerat=  sum(difmiu.*numerat1');
            nProb=(exp(numerat))/(denomin);
            [rNprob,cNprob]=size(nProb);
            nProb=reshape(nProb,rowX,colX,1);
            %if (rNprob~=1)&levX==1 nProb=diag(nProb); end
            return;
        end
    end
    if exist('difmiu')
        invS=inv(sigmax); detS=det(sigmax);
        denomin=sqrt(detS)*(2*pi)^(a/2);
        numerat=-0.5*(difmiu)'*invS*(difmiu);
        nProb=(exp(numerat))/(denomin);
        [rNprob,cNprob]=size(nProb);
        if levX==2 nProb=reshape(diag(nProb),rowX,colX,1); end
        if (rNprob~=1)&levX==1 nProb=diag(nProb); end
    else warning('case not covered ');  nProb=[];   end
end    



% % used in bayClass.m to calculate normal probability
% if iscell(x)
%     xLength=length(x);
%     x2=zeros(xLength,1);
%     for counter=1:xLength
%         x2(counter)=x{counter};
%     end
%     clear x;
%     x=x2;
%     clear x2;
% end
% if length(x)==32; x=x([1:22 24 25 28:32]); end
% if length(miux)==32 miux=miux([1:22 24 25 28:32]); end
% [a,b]=size(S);   if (a==32&b==32) S=S([1:22 24 25 28:32],[1:22 24 25 28:32]); a=a-3; end
% [c,d]=size(x);   if (c==1&d~=1) x=x'; end;
% [c,d]=size(miux); if (c==1&d~=1) miux=miux'; end;
% %invS=S;
% 
% difmiu=x-miux;
% %invS=1./(S);
% %invS=inv(S);
% %denomin=sqrt(det(S))*(2*pi)^(a/2);
% denomin=sqrt(detS)*(2*pi)^(a/2);
% 
% numerat=-0.5*(difmiu)'*invS*(difmiu);
% nProb=(exp(numerat))/(denomin);
% %nProb=numerat;
