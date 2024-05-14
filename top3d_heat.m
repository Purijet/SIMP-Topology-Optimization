clear; close all; clc
top3D(40,40,5,0.2,3.0,1.4)

% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
function top3D(nelx,nely,nelz,volfrac,penal,rmin)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
k0 = 1;           % Good thermal conductivity
kmin = 1e-3;      % Poor thermal conductivity
% USER-DEFINED FIXED DOFs
il=nelx/2-nelx/20:nelx/2+nelx/20; jl=nely; kl=0:nelz;                                   % nelx/2-nelx/20:nelx/2+nelx/20; jl=nely; kl=0:nelz
fixedxy = il*(nely+1)*(nely+1-jl);                                                      % il*(nely+1)*(nely+1-jl)
fixednid = repmat(fixedxy',size(kl)) + repmat(kl*(nelx+1)*(nely+1),size(fixedxy,2),1);
fixeddof = reshape(fixednid,[],1);
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = (nelx+1)*(nely+1)*(nelz+1);
F = sparse(1:ndof,1,-0.01,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(k0);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = nodeids(:)+1;
edofMat = repmat(edofVec,1,8)+ ...
    repmat([0 nely + [1 0] -1 ...
    (nely+1)*(nelx+1)+[0 nely + [1 0] -1]],nele,1);
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
x = repmat(volfrac,[nely,nelx,nelz]);
xPhys = x; 
loop = 0; 
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    % FE-ANALYSIS
    sK = reshape(KE(:)*(kmin+(1-kmin)*xPhys(:)'.^penal),8*8*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((kmin+(1-kmin)*xPhys.^penal).*ce)));
    dc = -penal*(k-kmin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx,nelz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);  
    dv(:) = H*(dv(:)./Hs);
    % OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        xPhys(:) = (H*xnew(:))./Hs;
        if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    % PLOT DENSITIES
    if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
end
clf; display_3D(xPhys);
end


% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(k)
A1 = 4*eye(2); A2 = -eye(2);
A3 = fliplr(A2); A4 = -ones(2);
KE1 = [A1 A2; A2 A1];
KE2 = [A3 A4; A4 A3];
KE = 1/12*k*[KE1 KE2; KE2 KE1];
end
% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
% =========================================================================
% === This code was written by K Liu and A Tovar, Dept. of Mechanical   ===
% === Engineering, Indiana University-Purdue University Indianapolis,   ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === The code is intended for educational purposes, and the details    ===
% === and extensions can be found in the paper:                         ===
% === K. Liu and A. Tovar, "An efficient 3D topology optimization code  ===
% === written in Matlab", Struct Multidisc Optim, 50(6): 1175-1196, 2014, =
% === doi:10.1007/s00158-014-1107-x                                     ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://www.top3dapp.com/             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, a ===