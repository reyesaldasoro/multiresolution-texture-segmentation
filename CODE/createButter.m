function butterf=createButter(rows,cols,levs,numOrient);
delta=1e-10;
numFeats=1;

if levs==1
    %----- 2D  Butterfly
    %     %----- 7 x 7   with 8 orientation
    butterf(:,:,:,1)=repmat([3 0 0 0 0 0 1;3 3 0 0 0 1 1;3 3 3 0 1 1 1;3 3 3 2 1 1 1;3 3 3 0 1 1 1;3 3 0 0 0 1 1;3 0 0 0 0 0 1]',[1 1 numFeats]);
    butterf(:,:,:,4)=repmat([3 3 3 3 3 0 0;3 3 3 3 0 0 1;3 3 3 0 0 1 1;3 3 3 2 1 1 1;3 3 0 0 1 1 1;3 0 0 1 1 1 1;0 0 1 1 1 1 1],[1 1 numFeats]);
    butterf(:,:,:,3)=repmat([3 3 3 3 0 0 0;3 3 3 3 0 0 0;3 3 3 3 0 0 0;3 3 3 2 1 1 1;0 0 0 1 1 1 1;0 0 0 1 1 1 1;0 0 0 1 1 1 1],[1 1 numFeats]);
    butterf(:,:,:,5)=repmat(butterf(:,:,1,1)',[1 1 numFeats]);
    butterf(:,:,:,2)=repmat(butterf(:,:,1,4)',[1 1 numFeats]);
    butterf(:,:,:,6)=butterf(:,end:-1:1,:,4);
    butterf(:,:,:,7)=butterf(:,end:-1:1,:,3);
    butterf(:,:,:,8)=butterf(:,end:-1:1,:,2);
    if numOrient==2
        butterf=butterf(:,:,:,[1 5]);
    elseif numOrient==4
        butterf=butterf(:,:,:,[1 3 5 7]);
    end
    if rows==5
        butterf=butterf(2:6,2:6,:,:);
    elseif rows==3
        butterf=butterf(3:5,3:5,:,:);
    end
else
    %----- 3D  Butterfly
    x=ones(rows,1)*linspace(-1,1,cols);y=x';
    X=prod3d(x,ones(levs,1));
    Y=prod3d(y,ones(levs,1));
    Z=prod3d(ones(rows,cols),linspace(-1,1,cols));
    R=sqrt(X.^2+Y.^2+Z.^2);R=R./(max(max(max(R))));
    thetaAxis=acos(Z./(R+delta));
    fiWidth=1.1*pi/4;
    fiAxis=-pi/2+(X<0)*pi+atan(Y./(X+delta));
    fiCentr=pi/2;
    q2=(((fiAxis>(fiCentr-fiWidth)).*(fiAxis<(fiCentr+fiWidth))));
    fiCentr=3*pi/4;thetaCentr=0;
    q3=(((fiAxis>(fiCentr-fiWidth)).*(fiAxis<(fiCentr+fiWidth))));
    q4=((permute(q2,[3 2 1])));
    q5=permute(q4,[2 1 3]);
    
    pAxis=q2&q4;
    pAxis3=pAxis(:,end:-1:1,:);
    pPlane=q3&q4&q5;
    pPlane3=pPlane(end:-1:1,end:-1:1,:);
    pQuadrant=(X>0)&(Y>0)&(Z>0);
    butterf=zeros(rows,cols,levs,numOrient,numOrient);
    %butterf(..x..,..y..,..z..,..FI..,..THETA..);
    %----- angles => dimension
    %----- 0        1       180     5
    %----- 45       2       215     6
    %----- 90       3       270     7
    %----- 135      4       305     8
    %----- Each butterfly pair will be symetric with the values 1 / 3
    %----- three along the axes  X - cols  Y - rows    Z - levs
    if (numOrient==2)
        % fi=90 theta=90  Y AXIS
        butterf(:,:,:,2,2)= (pAxis) + 3*(pAxis3) ;
        % fi=0 theta=90   X AXIS
        %butterf(:,:,:,1,3)= permute(pAxis,[2 1 3]) + 3*permute(pAxis3,[2 1 3])   ;
        butterf(:,:,:,1,2)= permute(butterf(:,:,:,2,2),[2 1 3]);
        % fi=0 theta=0    Z AXIS
        %butterf(:,:,:,1,1)= permute(pAxis,[3 2 1]) + 3*permute(pAxis3,[3 2 1]) ;
        butterf(:,:,:,1,1)= permute(butterf(:,:,:,2,2),[1 3 2]);
        %----- six along the planes
    end
    if numOrient==4
        % fi=90 theta=90  Y AXIS
        butterf(:,:,:,3,3)= (pAxis) + 3*(pAxis3) ;
        % fi=0 theta=90   X AXIS
        %butterf(:,:,:,1,3)= permute(pAxis,[2 1 3]) + 3*permute(pAxis3,[2 1 3])   ;
        butterf(:,:,:,1,3)= permute(butterf(:,:,:,3,3),[2 1 3]);
        % fi=0 theta=0    Z AXIS
        %butterf(:,:,:,1,1)= permute(pAxis,[3 2 1]) + 3*permute(pAxis3,[3 2 1]) ;
        butterf(:,:,:,1,1)= permute(butterf(:,:,:,3,3),[1 3 2]);
        butterf(:,:,:,2,1)=butterf(:,:,:,1,1);
        butterf(:,:,:,3,1)=butterf(:,:,:,1,1);
        butterf(:,:,:,4,1)=butterf(:,:,:,1,1);
        
        %----- six along the planes
        % fi=45 theta=90   XY 
        butterf(:,:,:,2,3)= (pPlane) + 3*(pPlane3) ;
        % fi=135 theta=90  XY 
        butterf(:,:,:,4,3)= permute(butterf(:,:,:,2,3),[3 2 1]) ;
        % fi=0 theta=45    XZ
        butterf(:,:,:,1,2)= permute(butterf(:,:,:,2,3),[1 3 2])   ;
        % fi=0 theta=135   XZ
        butterf(:,:,:,1,4)= butterf(:,:,end:-1:1,1,2) ;
        % fi=90 theta=45   YZ 
        butterf(:,:,:,3,2)= permute(butterf(:,:,:,2,3),[3 1 2])   ;
        % fi=90 theta=135  YZ 
        butterf(:,:,:,3,4)=  butterf(:,:,end:-1:1,1,4) ;
        %-----four along  ....the corners
        % fi=45 theta=45    1 1 1 
        butterf(:,:,:,2,2)= ((X>0)&(Y>0)&(Z>0)) +3*((X<0)&(Y<0)&(Z<0))  ;
        % fi=135 theta=45   -1 1 1 
        butterf(:,:,:,4,2)= ((X<0)&(Y>0)&(Z>0)) +3*((X>0)&(Y<0)&(Z<0))  ;
        % fi=215 theta=45   -1-1 1 
        butterf(:,:,:,2,4)= ((X<0)&(Y<0)&(Z>0)) +3*((X>0)&(Y>0)&(Z<0))  ;  % 6,2
        % fi=305 theta=45    1-1 1 
        butterf(:,:,:,4,4)= ((X>0)&(Y<0)&(Z>0)) +3*((X<0)&(Y>0)&(Z<0))  ;  % 8,2
    end
    if (numOrient~=2)&(numOrient~=4)
        warning('For 3D only 2 or 4 orientations is possible');
        butterf=[];
    end
end
    %----- butterflies
    % if numOrient==4
    %     % %----- 3 x 3    with 4 orientations
    %     butterf(:,:,1)=[1 0 1;1 1 1;1 0 1 ]; butterf(:,:,2)=[1 1 0;1 1 1;0 1 1 ];
    %     butterf(:,:,3)=[1 1 1;0 1 0;1 1 1 ]; butterf(:,:,4)=[0 1 1;1 1 1;1 1 0 ];
    % else
    %     %----- 5 x 5   with 8 orientation
    %     butterf(:,:,1)=[3 0 0 0 1;3 3 0 1 1;3 3 2 1 1;3 3 0 1 1;3 0 0 0 1]';
    %     butterf(:,:,4)=[3 3 3 0 0;3 3 0 0 1;3 3 2 1 1;3 0 0 1 1;0 0 1 1 1];
    %     butterf(:,:,3)=[3 3 3 0 0;3 3 3 0 0;3 3 2 1 1;0 0 1 1 1;0 0 1 1 1];
    %     butterf(:,:,2)=butterf(:,:,4)';butterf(:,:,5)=butterf(:,:,1)';
    %     butterf(:,:,6)=butterf(:,end:-1:1,4);butterf(:,:,7)=butterf(:,end:-1:1,3);
    %     butterf(:,:,8)=butterf(:,end:-1:1,2);
    % end
    %     if numOrient==4
    %         % %----- 3 x 3    with 4 orientations
    %         butterf(:,:,:,1)=repmat([3 3 3;0 2 0;1 1 1 ],[1 1 numFeats]); 
    %         butterf(:,:,:,2)=repmat([3 3 0;3 2 1;0 1 1 ],[1 1 numFeats]);
    %         butterf(:,:,:,3)=repmat([3 0 1;3 2 1;3 0 1 ],[1 1 numFeats]); 
    %         butterf(:,:,:,4)=repmat([0 3 3;1 2 3;1 1 0 ],[1 1 numFeats]);
    %     else
    %         %     %----- 5 x 5   with 8 orientation
    %         butterf(:,:,:,1)=repmat([3 0 0 0 1;3 3 0 1 1;3 3 2 1 1;3 3 0 1 1;3 0 0 0 1]',[1 1 numFeats]);
    %         butterf(:,:,:,4)=repmat([3 3 3 0 0;3 3 0 0 1;3 3 2 1 1;3 0 0 1 1;0 0 1 1 1],[1 1 numFeats]);
    %         butterf(:,:,:,3)=repmat([3 3 3 0 0;3 3 3 0 0;3 3 2 1 1;0 0 1 1 1;0 0 1 1 1],[1 1 numFeats]);
    %         %butterf(:,:,:,2)=butterf(:,:,:,4)';butterf(:,:,:,5)=butterf(:,:,:,1)';
    %         butterf(:,:,:,2)=repmat([3 3 3 0 0;3 3 0 0 1;3 3 2 1 1;3 0 0 1 1;0 0 1 1 1]',[1 1 numFeats]);
    %         butterf(:,:,:,5)=repmat([3 0 0 0 1;3 3 0 1 1;3 3 2 1 1;3 3 0 1 1;3 0 0 0 1],[1 1 numFeats]);
    %         butterf(:,:,:,6)=butterf(:,end:-1:1,:,4);butterf(:,:,:,7)=butterf(:,end:-1:1,:,3);
    %         butterf(:,:,:,8)=butterf(:,end:-1:1,:,2);
    %         %     %----- 7 x 7   with 8 orientation
    %         %     butterf(:,:,:,1)=repmat([3 0 0 0 0 0 1;3 3 0 0 0 1 1;3 3 3 0 1 1 1;3 3 3 2 1 1 1;3 3 3 0 1 1 1;3 3 0 0 0 1 1;3 0 0 0 0 0 1]',[1 1 numFeats]);
    %         %     butterf(:,:,:,4)=repmat([3 3 3 3 3 0 0;3 3 3 3 0 0 1;3 3 3 0 0 1 1;3 3 3 2 1 1 1;3 3 0 0 1 1 1;3 0 0 1 1 1 1;0 0 1 1 1 1 1],[1 1 numFeats]);
    %         %     butterf(:,:,:,3)=repmat([3 3 3 3 0 0 0;3 3 3 3 0 0 0;3 3 3 3 0 0 0;3 3 3 2 1 1 1;0 0 0 1 1 1 1;0 0 0 1 1 1 1;0 0 0 1 1 1 1],[1 1 numFeats]);
    %         %     butterf(:,:,:,5)=repmat(butterf(:,:,1,1)',[1 1 numFeats]);
    %         %     butterf(:,:,:,2)=repmat(butterf(:,:,1,4)',[1 1 numFeats]);
    %         %     butterf(:,:,:,6)=butterf(:,end:-1:1,:,4);
    %         %     butterf(:,:,:,7)=butterf(:,end:-1:1,:,3);
    %         %     butterf(:,:,:,8)=butterf(:,end:-1:1,:,2);
    %     end

