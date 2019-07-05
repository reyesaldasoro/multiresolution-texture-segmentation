function [redData]=reduceu(data,numReductions)
% function [redData]=reduceu(data)
%----------------------------------------------------------------
%------ REDUCEU  function to reduce data in uniform levels  -----
%------          receives an image and reduces it in quadtree ---
%------          reduction is ALWAYS in factor of 2 -------------
%----------------------------------------------------------------
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------        May  2001 ----------------------------
%----------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1    help reduceu;   redData=[];   return; end;
if ~exist('numReductions') numReductions=1; end
if ~(isa(data,'double')) data=double(data); end

if numReductions==0
    redData=data;
else
    if numReductions>1
        data=reduceu(data,numReductions-1);
    end
    
    [rows,cols,levels,otherDim]=size(data);
    
    if (levels>1)&(log2(levels)~=floor(log2(levels)))
             redData=convn(data,[1 1;1 1]);
            redData=redData(2:2:end,2:2:end,:,:)/4;
       
        return
    end
    %------ !!!!restriction!!! the data sizes a power of 2
%    if ((log2(rows)~=floor(log2(rows)))|(log2(cols)~=floor(log2(cols))))
        if levels==1
            %help reduceu;    warning('Data dimensions must be a power of 2');
            redData=conv2(data,[1 1;1 1]);
            redData=redData(2:2:end,2:2:end)/4;
        else
            redData=convn(data,ones(2,2,2));
            redData=redData(2:2:end,2:2:end,2:2:end)/8;

        end
        % redData=[];
        return;
%    end
    
    
%     %------ matrix to multiply in order to reduce the data ---------
%     %------ basic matrix, Identity modified   ----------------------
%     Imod=[1;1];                                                 %reduces 2 to 1
%     Imodn=setImod(Imod,rows);
%     
%     %-------------------general revision about the sizes ---------------------
%     if (rows<2)|(cols<2)
%         redData=(data'*Imodn/2)';
%     else
%         %---- determine if it is 4Dims
%         if otherDim>1
%             for counterDims=1:otherDim
%                 %determine 2D and 3D Cases
%                 if levels==1
%                     redData(:,:,:,counterDims)=Imodn'*data(:,:,:,counterDims)*Imodn/4;
%                 else
%                     for k=1:2:levels
%                         redData(:,:,(k+1)/2,counterDims)=(Imodn'*(data(:,:,k,counterDims)+data(:,:,k+1,counterDims))*Imodn)/8;
%                     end
%                 end
%             end
%             
%         else
%             %determine 2D and 3D Cases
%             if levels==1
%                 redData=Imodn'*data*Imodn/4;
%             else
%                 for k=1:2:levels
%                     redData(:,:,(k+1)/2)=(Imodn'*(data(:,:,k)+data(:,:,k+1))*Imodn)/8;
%                 end
%             end
%         end
%     end
end