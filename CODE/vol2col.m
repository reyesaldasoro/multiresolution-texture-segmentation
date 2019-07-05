function dataCol=vol2col(varargin)
%VOL2COL Rearrange volume blocks into columns in a sliding fashion.
%   B = IM2COL(A,[M N P]) converts each sliding M-by-N-by-P
%   block of A into a column of B, with no zero padding. B has
%   M*N*P rows and will contain as many columns as there are M-by-N-by-P
%   neighborhoods in A. If the size of A is [MM NN PP], then the
%   size of B is (M*N)-by-((MM-M+1)*(NN-N+1)*(PP-P+1)). Each column of B
%   contains the neighborhoods of A reshaped as NHOOD(:), where
%   NHOOD is a matrix containing an M-by-N-by-P neighborhood of
%   A. VOL2COL orders the columns of B so that they can be
%   reshaped to form a matrix in the normal way. For example,
%   suppose you use a function, such as SUM(B), that returns a
%   scalar for each column of B. You can directly store the
%   result in a matrix of size (MM-M+1)-by-(NN-N+1) using these
%   calls: 
%
%        B = im2col(A,[M N],'sliding');
%        C = reshape(sum(B),MM-M+1,NN-N+1);
%----------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro based on IM2COL
%------             PHD     the University of Warwick
%------  Supervisor :   Abhir Bhalerao    -----------
%------  19 September 2002 --------------------------
%----------------------------------------------------

%------ no input data is received, error -------------------------
if nargin<1   help vol2col;    b=[];   return; end;

data=varargin{1};              %----- Volume to be transformed
block=varargin{2};          %----- Block size of the transformation
if nargin>2
   boundary=varargin{3};   %----- Boundary for a reduced version
end

[ma,na,pa] = size(data);

if pa==1 
    if nargin==2
        if length(block)==3 block=block(1:2); end
        dataCol=im2col(data,block);      %----- In case is 2D use IM2COL
    else
         dataCol=im2colRed(data,block,boundary);      %----- In case is 2D use IM2COL
    end
else
    m = block(1);
    if length(block)==1         n=1;p=1;
    elseif length(block)==2     n = block(2); p=1;
    else                        n = block(2); p=block(3);
    end
	if nargin>2
        boundary=varargin{3};   %----- Boundary for a reduced version
        if all(size(boundary)==size(data))
            boundary2=boundary(1+floor(m/2):end-floor(m/2),1+floor(n/2):end-floor(n/2),1+floor(p/2):end-floor(p/2));
        else
            boundary2=boundary;
        end
	end
        
    if any([ma na pa] < [m n p]) % if neighborhood is larger than image
        dataCol = []; %zeros(m*n,0);
        return
    end
    
   
    % Create Hankel-like indexing sub matrix.
    mc = block(1); nc = ma-m+1; nn = na-n+1;        %nc and nn get the number of effective pixels to use 
    cidx = (0:mc-1)'; ridx = 1:nc;                  %cidx column vector ridx row vector
    %% cidx(:,ones(nc,1)) is column form increase   [0 0 ..;1 1..; 2 2.. ]
    %% ridx(:,ones(nc,1)) is a row form increase    [0 1 2 ..; 0 1 2..; 0 1 2.. ]
    t = cidx(:,ones(nc,1)) + ridx(ones(mc,1),:);    % Hankel Subscripts
    %% the sum is a combination of the indexes  t = [1 2 3 ..; 2 3 4.. ; 3 4 5..]
    colIndex = zeros(mc*n,nc);                            %% colIndex [ numNeighbours , numPixelsInOneColumn]
    rows = [1:mc];
    for i=0:n-1,
        colIndex(i*mc+rows,:) = t+ma*i;
    end
    rowColIndex = zeros(mc*n,nc*nn);
    cols = 1:nc;
    %%% rowColIndex becomes an index to extract the values from "a" it is a linear
    %%% (1D) index that looks for the values in the matrix
    %%% the cycle adds the next columns (colIndex holds the first one)

    
    for j=0:nn-1,
        rowColIndex(:,j*nc+cols) = colIndex+ma*j;
    end
    pc=block(3);
    pp=pa-p+1;
    rows2=1:size(rowColIndex,1);
    rowColLevIndex=zeros(mc*n*p,nc*nn);
    for i=0:pc-1
        rowColLevIndex(i*m*n+rows2,:)=rowColIndex+ma*na*i;
    end
    %rowColLevTotIndex=zeros(mc*n*p,nc*nn*pp);
    rowColLevTotIndex=[];
    cols2=1:size(rowColLevIndex,2);
    for j=0:pp-1    
%        rowColLevTotIndex(:,j*pp*nc +cols2)=rowColLevIndex+j*ma*na;
            %----- reduction is achieved by not including all the columns
            %----- of the penultimate plane
            if exist('boundary')
                bounPlane=boundary2(:,:,j+1);
                reducedPlane=rowColLevIndex(:,bounPlane(:)~=0);
            else
                reducedPlane=rowColLevIndex;
            end
            rowColLevTotIndex=[rowColLevTotIndex reducedPlane+j*ma*na];
            %rowColLevTotIndex(:,j*nn*nc +cols2)=rowColLevIndex+j*ma*na;
    end
    
    dataCol = data(rowColLevTotIndex);
end