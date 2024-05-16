clear; close all; clc;
addpath('./Mesh_voxelisation')

% Read The CAD Model
[bcc.facets, bcc.vertices, bcc.normals] = stlread2('BCCunit.stl');
figure
trisurf(bcc.facets,bcc.vertices(:,1),bcc.vertices(:,2),bcc.vertices(:,3),'FaceColor','black','FaceAlpha',0.1);
axis equal; axis tight; axis on; box on; view([30,30]);
CopyLattice(bcc.facets,bcc.vertices,3,3,3)

% FCCvert = [];
% [faceCell, vertCell, normCell] = stlread2('BCC.stl');
% figure
% scatter3(vertCell(:,1),vertCell(:,2),vertCell(:,3),'.')
% axis equal; axis tight; axis on; box on; view([30,30]);
% figure
% trimesh(faceCell,vertCell(:,1),vertCell(:,2),vertCell(:,3))
% axis equal; axis tight; axis on; box on; view([30,30]);

% Read The CAD Model
[stl.facets, stl.vertices, stl.normals] = stlread2('Femoral component.stl');
FindOverhang(stl.facets, stl.vertices, stl.normals, 0, 0)
[SupportGrid, PassGrid] = VoxelizeSTL(stl.vertices,0,0);
[nelx, nely, nelz] = size(SupportGrid);

VertX = []; VertY = []; VertZ = []; Face = []; num = 0;
for i = 1:nelx
    for j = 1:nely
        for k = 1:nelz
            if SupportGrid(i,j,k) == 1
                VertX = cat(1,VertX,bcc.vertices(:,1)/(max(bcc.vertices(:,1))-min(bcc.vertices(:,1))) + i);
                VertY = cat(1,VertY,bcc.vertices(:,2)/(max(bcc.vertices(:,2))-min(bcc.vertices(:,2))) + j);
                VertZ = cat(1,VertZ,bcc.vertices(:,3)/(max(bcc.vertices(:,3))-min(bcc.vertices(:,3))) + k);
                Face = cat(1,Face,bcc.facets + num*length(bcc.vertices));
                num = num + 1;
            end
        end
    end
end
Vert = [VertX VertY VertZ];
figure
scatter3(Vert(:,1),Vert(:,2),Vert(:,3),'.')
axis equal; axis tight; axis on; box on; view([30,30]);
figure
trimesh(Face,Vert(:,1),Vert(:,2),Vert(:,3))
axis equal; axis tight; axis on; box on; view([30,30]);

% Plot overhang area of the part
function FindOverhang(Facets, Vertices, Normals, alpha, beta)
    Rox = [1 0 0 0; 0 cos(alpha) sin(alpha) 0; 0 -sin(alpha) cos(alpha) 0; 0 0 0 1];
    Roy = [cos(beta) 0 -sin(beta) 0; 0 1 0 0; sin(beta) 0 cos(beta) 0; 0 0 0 1];
    Vertices(:,end+1) = ones(1,length(Vertices))'; Vertices = Vertices*Rox*Roy;
    T = Facets; x = Vertices(:,1); y = Vertices(:,2); z = Vertices(:,3);
    Normals(:,end+1) = ones(1,length(Normals))'; F = Normals*Rox*Roy;
    T_hang = []; T_nohang = []; Area = [];
    for i = 1:size(T,1)
        if F(i,3) < -0.5
            T_hang(end+1,:) = T(i,:);
            Area(end+1) = triangle(Vertices(T(i,1),1:3),Vertices(T(i,2),1:3),Vertices(T(i,3),1:3));
        else
            T_nohang(end+1,:) = T(i,:);
        end
    end
    figure
    trisurf(T_nohang,x,y,z,'FaceColor','black','FaceAlpha',0.1);
    axis equal;  axis tight; axis on; box on; view([30,30]);
    hold on
    trisurf(T_hang,x,y,z,'FaceColor','red','FaceAlpha',0.8);
    fprintf('The areas in order are %.2f\n', sum(Area));
end

% Heron Formula
function [area] =  triangle(A,B,C)
    a = sqrt((A(1)-B(1))^2+(A(2)-B(2))^2+(A(3)-B(3))^2);
    b = sqrt((B(1)-C(1))^2+(B(2)-C(2))^2+(B(3)-C(3))^2);
    c = sqrt((C(1)-A(1))^2+(C(2)-A(2))^2+(C(3)-A(3))^2);
    s = (a+b+c)/2;
    area = sqrt(s*(s-a)*(s-b)*(s-c));
end

function CopyLattice(face,vert,numX,numY,numZ)
    VertX = []; VertY = []; VertZ = []; Face = []; num = 0;
    for i = 1:numX
        for j = 1:numY
            for k = 1:numZ
                VertX = cat(1,VertX,vert(:,1)/(max(vert(:,1))-min(vert(:,1))) + i);
                VertY = cat(1,VertY,vert(:,2)/(max(vert(:,2))-min(vert(:,2))) + j);
                VertZ = cat(1,VertZ,vert(:,3)/(max(vert(:,3))-min(vert(:,3))) + k);
                Face = cat(1,Face,face + num*length(vert));
                num = num + 1;
            end
        end
    end
    Vert = [VertX VertY VertZ];
    figure
    scatter3(Vert(:,1),Vert(:,2),Vert(:,3),'.')
    axis equal; axis tight; axis on; box on; view([30,30]);
    figure
    trimesh(Face,Vert(:,1),Vert(:,2),Vert(:,3))
    axis equal; axis tight; axis on; box on; view([30,30]);
end

% Voxelize the part and generate voxel support structures
function [SupportGrid, PassGrid] = VoxelizeSTL(Vertices, alpha, beta)
    Rox = [1 0 0 0; 0 cos(alpha) sin(alpha) 0; 0 -sin(alpha) cos(alpha) 0; 0 0 0 1];
    Roy = [cos(beta) 0 -sin(beta) 0; 0 1 0 0; sin(beta) 0 cos(beta) 0; 0 0 0 1];
    V = Vertices; V(:,end+1) = ones(1,length(Vertices))'; V = V*Rox*Roy;
    voxel_size = 3;
    xco = zeros(3,size(V,1)/3);
    yco = zeros(3,size(V,1)/3);
    zco = zeros(3,size(V,1)/3);
    for i=1:size(V)/3
        xco(1,i)=V(3*i-2,1);
        xco(2,i)=V(3*i-1,1);
        xco(3,i)=V(3*i,1);
    
        yco(1,i)=V(3*i-2,2);
        yco(2,i)=V(3*i-1,2);
        yco(3,i)=V(3*i,2);
    
        zco(1,i)=V(3*i-2,3);
        zco(2,i)=V(3*i-1,3);
        zco(3,i)=V(3*i,3);
    end
    n_grid_x = round(abs((max(V(:,1))-min(V(:,1)))/voxel_size));
    n_grid_y = round(abs((max(V(:,2))-min(V(:,2)))/voxel_size));
    n_grid_z = round(abs((max(V(:,3))-min(V(:,3)))/voxel_size));
    [OUTPUTgrid,x,y,z] = VOXELISE(n_grid_x,n_grid_y,n_grid_z,xco,yco,zco);
    OutputGrid = zeros(n_grid_x,n_grid_y,n_grid_z+5); OutputGrid(:,:,6:n_grid_z+5)=OUTPUTgrid;
    SupportGrid = zeros(n_grid_x,n_grid_y,n_grid_z+5);
    PassGrid = ones(n_grid_x,n_grid_y,n_grid_z+5);
    figure
    hold on;
    axis equal; view([30,30]);
    Support = 0;
    face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
    hx=x(2)-x(1);
    hy=y(2)-y(1);
    hz=z(2)-z(1);
    for i=1:n_grid_x
        for j=1:n_grid_y
            partVoxelFound = false;
            for k=flip(1:n_grid_z)
                if (OutputGrid(i,j,k)==1)
                    partVoxelFound = true;
                    vert=[x(i) y(j) z(k); x(i) y(j)-hy z(k); x(i)+hx y(j)-hy z(k);...
                        x(i)+hx y(j) z(k); x(i) y(j) z(k)+hz; x(i) y(j)-hy z(k)+hz;...
                        x(i)+hx y(j)-hy z(k)+hz; x(i)+hx y(j) z(k)+hz];
                    patch('Faces',face,'Vertices',vert,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.3);
                    hold on;
                end
                if (OutputGrid(i,j,k)==0 && partVoxelFound==true)
                    SupportGrid(i,j,k) = 1; PassGrid(i,j,k) = 0;
                    Support = Support + voxel_size^3;
                    vert=[x(i) y(j) z(k); x(i) y(j)-hy z(k); x(i)+hx y(j)-hy z(k);...
                        x(i)+hx y(j) z(k); x(i) y(j) z(k)+hz; x(i) y(j)-hy z(k)+hz;...
                        x(i)+hx y(j)-hy z(k)+hz; x(i)+hx y(j) z(k)+hz];
                    patch('Faces',face,'Vertices',vert,'FaceColor',[0.7 0 0],'FaceAlpha',0.3);
                    hold on;
                end
            end
        end
    end
    fprintf('The support volume in order are %.2f\n', Support);
end