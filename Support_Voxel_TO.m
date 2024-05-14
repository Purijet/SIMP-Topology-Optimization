clear; close all; clc;
addpath('./Mesh_voxelisation')

% Read The CAD Model
[stl.facets, stl.vertices, stl.normals] = stlread2('Femoral component.stl');
FindOverhang(stl.facets, stl.vertices, stl.normals, 0, 0)
[SupportGrid, PassGrid] = VoxelizeSTL(stl.vertices,0,0);
[nelx, nely, nelz] = size(SupportGrid);
for i = 1:nelx
    for j = 1:nely
        for k = 1:4
            if (mod(i,2)==1 || mod(j,2)==1)
                SupportGrid(i,j,k) = 0;
                PassGrid(i,j,k) = 1;
            end
        end
    end
end
figure
display_3D(PassGrid)
% Top3dSTL_v3('VoxelSupport.STL',SupportGrid,'Format','binary')

% Rotate 30 degree along X-axis and 0 degree along Y-axis
% FindOverhang(stl.facets, stl.vertices, stl.normals, pi/6, 0)
% VoxelizeSTL(stl.vertices,pi/6,0)

% Rotate 50 degree along X-axis and -30 degree along Y-axis
% FindOverhang(stl.facets, stl.vertices, stl.normals, pi*5/18, -pi/6)
% VoxelizeSTL(stl.vertices, pi*5/18, -pi/6)

fixedxyz = [];
for i = 1:nelx
    for j = 1:nely
        for k = 1:nelz
            if (SupportGrid(i,j,k)==1)  % Overhang Detection
                fixedxyz(end+1,:) = [i,j,k]; % il(end+1) = i; jl(end+1) = j; kl(end+1) = k;
            end
        end
    end
end
fixed_bc = []; 
for i = 1:nelx
    for j = 1:nely
        load_bc_find = [];
        for e = 1:length(fixedxyz)
            if (fixedxyz(e,1) == i && fixedxyz(e,2) == j)
                load_bc_find(end+1,:) = fixedxyz(e,:);
            end
        end
        if isempty(load_bc_find),continue; else,fixed_bc(end+1,:) = load_bc_find(end,:); end
    end
end
figure
display_3D(SupportGrid);
hold on
scatter3(fixed_bc(:,1),fixed_bc(:,2),fixed_bc(:,3),'.','red');
axis equal; axis tight; box on; view([30,30]);

TOP3d_Heat(PassGrid,fixed_bc,nelx,nely,nelz);

%% === TOP3d Heat ===
function TOP3d_Heat(PassGrid,fixed_bc,nelx,nely,nelz)
    volfrac=0.3; penal=3.0; rmin=1.4;
    volfrac=volfrac*length(find(PassGrid(:)))/length(PassGrid(:));
    % USER-DEFINED LOOP PARAMETERS
    maxloop = 200;    % Maximum number of iterations
    tolx = 0.01;      % Terminarion criterion
    displayflag = 0;  % Display structure flag
    % USER-DEFINED MATERIAL PROPERTIES
    k0 = 1;           % Good thermal conductivity
    kmin = 1e-3;      % Poor thermal conductivity
    % USER-DEFINED LOAD DOFs
    il=fixed_bc(:,1); jl=fixed_bc(:,2); kl=fixed_bc(:,3);                                   % nelx/2-nelx/20:nelx/2+nelx/20; jl=nely; kl=0:nelz
    fixedxy = il*(nely+1).*(nely+1-jl);                                                      % il*(nely+1)*(nely+1-jl)
    fixednid = repmat(fixedxy',size(kl)) + repmat(kl*(nelx+1)*(nely+1),size(fixedxy,2),1);
    fixeddof = reshape(fixednid,[],1);
    % % PREPARE FINITE ELEMENT ANALYSIS
    nele = nelx*nely*nelz;
    ndof = (nelx+1)*(nely+1)*(nelz+1);
    alldof = reshape(1:ndof,[nely+1,nelx+1,nelz+1]);
    F = sparse(1:ndof,1,-0.01,ndof,1);
    U = zeros(ndof,1);
    freedofs = setdiff(1:ndof,fixeddof);
    KE = lk_H8(k0);
    nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
    nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
    nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
    nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
    % nodeids = reshape(1:(nely+1)*(nelx+1)*(nelz+1),nely+1,nelx+1,nelz+1);
    edofVec = nodeids(:)+1;
    edofMat = repmat(edofVec,1,8)+ repmat([0 nely + [1 0] -1 (nely+1)*(nelx+1)+[0 nely + [1 0] -1]],nele,1);
    iK = reshape(kron(edofMat,ones(8,1))',8*8*nele,1);
    jK = reshape(kron(edofMat,ones(1,8))',8*8*nele,1);
    % PREPARE FILTER
    iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for k1 = 1:nelz
        for i1 = 1:nelx
            for j1 = 1:nely
                e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
                for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                            e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                        end
                    end
                end
            end
        end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);
    % INITIALIZE ITERATION
    x = repmat(volfrac,[nelx,nely,nelz]);
    % row = []; col = []; axi = [];
    % for ely = 1:nely  % passive
    %     for elx = 1:nelx
    %         if sqrt((ely-nely/2.)^2+(elx-nelx/3.)^2) < nely/3.
    %             passive(ely,elx) = 1;
    %             row(end+1) = elx;
    %             col(end+1) = ely;
    %             axi(end+1) = 0;
    %         else
    %             passive(ely,elx) = 0;
    %         end
    %     end
    % end
    % passive = repmat(passive,[1,1,nelz]);
    x(find(PassGrid)) = 0;
    xPhys = x;
    disp([nelx,nely,nelz])
    loop = 0; 
    change = 1;
    loss_it = []; loss_obj = []; loss_vol = []; loss_ch = [];
    % START ITERATION
    while change > tolx && loop < maxloop
        loop = loop+1;
        % FE-ANALYSIS
        sK = reshape(KE(:)*(kmin+(1-kmin)*xPhys(:)'.^penal),8*8*nele,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nelx,nely,nelz]);
        c = sum(sum(sum((kmin+(1-kmin)*xPhys.^penal).*ce)));
        dc = -penal*(k-kmin)*xPhys.^(penal-1).*ce;
        dv = ones(nelx,nely,nelz);
        % FILTERING AND MODIFICATION OF SENSITIVITIES
        dc(:) = H*(dc(:)./Hs);  
        dv(:) = H*(dv(:)./Hs);
        % OPTIMALITY CRITERIA UPDATE
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            % xnew(find(passive)) = 0;  % passive
            xPhys(:) = (H*xnew(:))./Hs;
            if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else, l2 = lmid; end
        end
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        % PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
        loss_it(end+1)=loop; loss_obj(end+1)=c; loss_vol(end+1)=mean(xPhys(:)); loss_ch(end+1)=change;
        % PLOT DENSITIES
        if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
    end
    figure
    display_3D(xPhys);
    
    figure
    plot(loss_it,loss_obj)
    title('Compliance')
    xlabel('Iteration')
    ylabel('Compliance (J)')
    figure
    plot(loss_it,loss_vol)
    title('Volume fraction')
    xlabel('Iteration')
    ylabel('Volume fraction')
    figure
    plot(loss_it,loss_ch)
    title('Loss change')
    xlabel('Iteration')
    ylabel('Loss change')
end

% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(k)
A1 = 4*eye(2); A2 = -eye(2);
A3 = fliplr(A2); A4 = -ones(2);
KE1 = [A1 A2; A2 A1];
KE2 = [A3 A4; A4 A3];
KE = 1/12*k*[KE1 KE2; KE2 KE1];
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

%% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nelx,nely,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = (j-1)*hy;
            if (rho(i,j,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y+hy z; x+hx y+hy z; x+hx y z; x y z+hz;x y+hy z+hz; x+hx y+hy z+hz;x+hx y z+hz];
                vert(:,[2 3]) = vert(:,[2 3]); vert(:,2,:) = vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(i,j,k)),0.2+0.8*(1-rho(i,j,k)),0.2+0.8*(1-rho(i,j,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end