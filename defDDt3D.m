function [D,Dt] = defDDt3D

D = @(U) ForwardD(U);
Dt = @(X,Y,Z) Dive(X,Y,Z);

function [Dux,Duy,Duz] = ForwardD(U)
% [ux,uy,uz] = D u

Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
Duz = [diff(U,1,3); U(:,:,1) - U(:,:,end)];


function DtXYZ = Dive(X,Y,Z)
% DtXYZ = D_1' X + D_2' Y + D_3' Z

DtXYZ = [X(:,end,:) - X(:,1,:), -diff(X,1,2)];
DtXYZ = DtXYZ + [Y(end,:,:) - Y(1,:,:); -diff(Y,1,1)];
DtXYZ = DtXYZ + [Z(:,:,end) - Z(:,:,1); -diff(Z,1,3)];
DtXYZ = DtXYZ(:);