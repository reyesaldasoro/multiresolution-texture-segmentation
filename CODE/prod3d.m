function XYZ=prod3d(x,y,z)
%function xyz=prod3d(x,y,z)

if exist('x') [rx,cx,lx]=size(x); else rx=1;cx=1;lx=1; end
if exist('y') [ry,cy,ly]=size(y); else ry=1;cy=1;ly=1; end
if exist('z') [rz,cz,lz]=size(z); else rz=1;cz=1;lz=1; end





if nargin==2; 
    XY=x; z=y; [rz,cz,lz]=size(z);
else
    x=reshape(x,max(rx,cx),1);y=reshape(y,1,max(ry,cy));
    XY=x*y;
end

for counZ=1:max(rz,cz)
    XYZ(:,:,counZ)=XY*(z(counZ));
end
