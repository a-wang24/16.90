% This script will refine an existing mesh for the
% cylinder flow problem using a list of triangles to refine.
% This list must be in the array RefineList and is simply
% the list of integer indices of the triangles to be refined.

% Check to see if RefineList has been defined
if (~exist('RefineList')),
  fprintf('ERROR: RefineList does not exist.\n');
  fprintf('       You must create a list of cells to refine.\n');
  return;
end

% Load in the cylinder geometry description
% using a Matlab Decomposed Geometry Matrix.
%
load cylgeom; % The Decomposed Geom. Matrix is in the
              % variable cylgeom after this load.

% Calculate nodal data by averaging cell data
Unode = zeros(Nv,4);
Anode = zeros(Nv,1);
for i = 1:Nt,
  for j = 1:3,
    Unode(tri2nod(j,i),:) = Unode(tri2nod(j,i),:) + A(i)*U(:,i)';
    Anode(tri2nod(j,i))   =   Anode(tri2nod(j,i)) + A(i);
  end
end
Unode(:,1) = Unode(:,1)./Anode;
Unode(:,2) = Unode(:,2)./Anode;
Unode(:,3) = Unode(:,3)./Anode;
Unode(:,4) = Unode(:,4)./Anode;

% Perform refinement
[cyl_p,cyl_e,cyl_t,Unode_new] = refinemesh(cylgeom,cyl_p,cyl_e,cyl_t,Unode,RefineList);

% Set up all of the meshing info
SetupMesh;

% Move new nodal data back to cells by averaging
U = zeros(4,Nt);
for i = 1:Nt,
  for j = 1:3,
    U(:,i)= U(:,i) + Unode_new(tri2nod(j,i),:)';
  end
end
U = (1/3)*U;

% Plot adapted mesh
figure(3);
pdeplot(cyl_p, cyl_e, cyl_t);
axis('equal');
