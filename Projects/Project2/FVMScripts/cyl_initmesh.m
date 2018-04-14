% Load in the cylinder geometry description
% using a Matlab Decomposed Geometry Matrix.
%
load cylgeom; % The Decomposed Geom. Matrix is in the
              % variable cylgeom after this load.

% Generate first mesh
[cyl_p,cyl_e,cyl_t] = initmesh(cylgeom);

% This first mesh is a little too coarse so generate
% a refinement.
[cyl_p,cyl_e,cyl_t] = refinemesh(cylgeom,cyl_p,cyl_e,cyl_t);

% Set up all of the meshing info
SetupMesh;

