function dataCol=im2colRed(data,block,boundary)
%function dataCol=im2colRed(data,block,boundary)
%IM2COLRED Rearrange image blocks into columns with only some of the pixels
%of the image modification by Constantino Carlos Reyes Aldasoro 13 Oct 2003
%   restrictions:   only sliding no distinct
%   B = IM2COL(A,[M N]) rearranges each ****sliding ******  O N L Y 
%   M-by-N block in the image A into a column of B. IM2COL pads A
%   with zeros, if necessary, so its size is an integer multiple
%   of M-by-N. If A = [A11 A12; A21 A22], where each Aij is
%   M-by-N, then B = [A11(:) A12(:) A21(:) A22(:)].
%
%   B = IM2COL(A,[M N],'sliding') converts each sliding M-by-N
%   block of A into a column of B, with no zero padding. B has
%   M*N rows and will contain as many columns as there are M-by-N
%   neighborhoods in A. If the size of A is [MM NN], then the
%   size of B is (M*N)-by-((MM-M+1)*(NN-N+1). Each column of B
%   contains the neighborhoods of A reshaped as NHOOD(:), where
%   NHOOD is a matrix containing an M-by-N neighborhood of
%   A. IM2COL orders the columns of B so that they can be
%   reshaped to form a matrix in the normal way. For example,
%   suppose you use a function, such as SUM(B), that returns a
%   scalar for each column of B. You can directly store the
%   result in a matrix of size (MM-M+1)-by-(NN-N+1) using these
%   calls: 
%
%        B = im2col(A,[M N],'sliding');
%        C = reshape(sum(B),MM-M+1,NN-N+1);
%
%   B = IM2COL(A,[M N]) uses the default block type of
%   'sliding'.
%
%   B = IM2COL(A,'indexed',...) processes A as an indexed image,
%   padding with zeros if the class of A is uint8 or uint16, or 
%   ones if the class of A is double.
%
%   Class Support
%   -------------
%   The input image A can be numeric or logical. The output matrix 
%   B is of the same class as the input image.
%
%   See also BLKPROC, COL2IM, COLFILT, NLFILTER.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.22 $  $Date: 2002/03/29 15:55:08 $

    [ma,na] = size(data);
    m = block(1); n = block(2);
    if any([ma na] < [m n]) % if neighborhood is larger than image
       dataCol = zeros(m*n,0);
       return
    end
    % Create Hankel-like indexing sub matrix.
    mc = block(1); nc = ma-m+1; nn = na-n+1;        %nc and nn get the number of effective pixels to use 
    cidx = (0:mc-1)'; ridx = 1:nc;                  %cidx column vector ridx row vector
    %% cidx(:,ones(nc,1)) is column form increase   [0 0 ..;1 1..; 2 2.. ]
    %% ridx(:,ones(nc,1)) is a row form increase    [0 1 2 ..; 0 1 2..; 0 1 2.. ]
    t = cidx(:,ones(nc,1)) + ridx(ones(mc,1),:);    % Hankel Subscripts
    %% the sum is a combination of the indexes  t = [1 2 3 ..; 2 3 4.. ; 3 4 5..]
    %% colIndex is keeping a whole column of the original data and it will be
    %% used to generate the column values (original data) into dataCol next columns 
    colIndex = zeros(mc*n,nc);                            %% colIndex [ numNeighbours , numPixelsInOneColumn]
    rows = [1:mc];
    for i=0:n-1,
        colIndex(i*mc+rows,:) = t+ma*i;
    end
    rowColIndex=[];%rowColIndex = zeros(mc*n,nc*nn);
    cols = 1:nc;
    %%% rowColIndex becomes an index to extract the values from "a" it is a linear
    %%% (1D) index that looks for the values in the matrix
    %%% the cycle adds the next columns (colIndex holds the first one)
    

    if all(size(boundary)==size(data))
        boundary2=boundary(1+floor(m/2):end-floor(n/2),1+floor(m/2):end-floor(n/2));
    else
        boundary2=boundary;
    end
    
    for j=0:nn-1
        reducedColumn=colIndex(:,boundary2(:,j+1)~=0);
        rowColIndex=[rowColIndex reducedColumn+ma*j]; %rowColIndex(:,j*nc+cols) = colIndexRed+ma*j;
    end
    dataCol = data(rowColIndex);

%  [ma,na] = size(a);   
% b=[];
% rSize=floor(block(1)/2);
% cSize=floor(block(2)/2);
% %counter determine the centre of the block
% for rCoun=rSize+1:ma-rSize
%     for cCoun=cSize+1:na-cSize
%         if boundary(rCoun-rSize,cCoun-cSize)~=0
%             cDims=cCoun-cSize:cCoun+cSize;
%             rDims=rCoun-rSize:rCoun+rSize;
%             x=a(rDims,cDims);
%             b=[b x(:)];
%         end
%     end
% end

    
