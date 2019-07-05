function [bhattm,b1,b2]=bhattaM(mi,mj,si,sj)
% function  bhattm]=bhattam_r(mi,mj,si,sj)
%------  function that calculates the Bhattacharyya distance measure as
%------  presented by coleman et al
%------  input :    mean i, j
%------             sigma i, j   ~= 0  
%------             sigma==0  ---> infinite distance i.e. immediate 
%------             differenciation and no need  will be assigned as -1
%------  output :   bhattacharyya distance
%-----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  8 August 2002 ----------------------------
%----------------------------------------------------

%-----------------------------------------------------------------
%------ no input data is received, error -------------------------
%-----------------------------------------------------------------
if nargin < 1  help bhattaM; bhattaM=[]; return; end;

delta1=1e-10;

 %
%-----------------------------------------------------------------
%------ get basic parameters of the data and revise dimensions
%-----------------------------------------------------------------
 [rowsmi,colsmi,levelsmi]=size(mi);        %----- dimensions mean i
 [rowsmj,colsmj,levelsmj]=size(mj);        %----- dimensions mean j

%-----------------------------------------------------------------
%------ input data is the samples, calculate the means and stds --
%-----------------------------------------------------------------
if (nargin==2)&~all(size(mj)==1)
    m_i=sum(sum(sum(mi)))/rowsmi/colsmi/levelsmi;   %----- mean first set
    m_j=sum(sum(sum(mj)))/rowsmj/colsmj/levelsmj;   %----- mean second set
    si=std2(mi);  %std(std(std(mi)));               *----- std first set
    sj=std2(mj);  %std(std(std(mj)));               %----- std second set
    mi=m_i; mj=m_j;
    clear m_i m_j;
end

if (all(size(mi)==1))&(all(size(mj)==1))    %------ this covers the 1D case
    %-----------------------------------------------------------------
    %------in case there are zeros in stds ---------------------------
    %-----------------------------------------------------------------
    if any(any(any(si==0)))
        tempsi=si;
        si(si==0)=1;
    end
    if any(any(any(sj==0)))
        tempsj=sj;
        sj(sj==0)=1;
    end
    %-----------------------------------------------------------------
    %------ bhattacharyya distance Measure ---------------------------
    %-----------------------------------------------------------------
    si2=si.^2;
    sj2=sj.^2;
    b1=0.25*log(((si2./sj2)+(sj2./si2)+2)/4);
    b2=0.25*(((mi-mj).^2)./(si2+sj2));
    bhattm=b1+b2;
    
    if exist('tempsi')    
        bhattm(tempsi==0)=-1;
    end
    if exist('tempsj')    
        bhattm(tempsj==0)=-1;
    end
else                                        %------ this covers the N-D case with covariance matrix!
    %----------------------------------------------------------------------------------
    %------ there are two cases, (a) two vectors of same dims, means and point --------
    %------ (b) one vector of means and a matrix 3D of images -------------------------
    %----------------------------------------------------------------------------------
    if (rowsmi==1|colsmi==1)
        [rowssi,colssi,levelssi]=size(si);        %----- dimensions mean i
        [rowssj,colssj,levelssj]=size(sj);        %----- dimensions mean j
        si2=si.^2;
        sj2=sj.^2;
        b1=0.25*log(((si2./(sj2+delta1))+(sj2./(si2+delta1))+2)/4);
        b2=0.25*(((mi-mj).^2)./(si2+sj2+delta1));
        bhattm=b1+b2;
    
    else
        
        [rowssi,colssi,levelssi]=size(si);        %----- dimensions mean i
        [rowssj,colssj,levelssj]=size(sj);        %----- dimensions mean j
        if ~(rowssi==rowssj&colssi==colssj&levelssi==levelssj) poolCovMatrix=si; 
        else     poolCovMatrix   = (si*si'+sj*sj')/2;
        end
        invCovMatrix    = inv(poolCovMatrix);
        
        if (rowsmi==rowsmj&colsmi==colsmj&levelsmi==levelsmj) 
            diffMius        = (mi-mj);
            b1=(1/8)*diffMius'*(invCovMatrix)*diffMius;
        else 
            diffMius=(reshape(mj,rowsmj*colsmj,rowsmi,1))-repmat(mi',rowsmj*colsmj,1);
            b1=(1/8)*sum((diffMius*invCovMatrix)'.*diffMius');
            b1=reshape(b1,rowsmj,colsmj,1);
        end
        b2=0.5*(log( det(poolCovMatrix) / (  sqrt(det(si)*det(sj))   ) ));
        bhattm=b1+b2;
    end
end