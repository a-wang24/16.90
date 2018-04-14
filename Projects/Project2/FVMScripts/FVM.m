% This Matlab script solves the two-dimensional Euler
% equations using a finite volume algorithm for the flow
% around a cylinder.
%

% If FVMrestart does not exist, set it to zero (assume not restarting)
if (~exist('FVMrestart')),
  FVMrestart = 0;
end

if (FVMrestart == 0),
  % Not a restart, need to initialize a bunch of data
  cyl_initmesh;

  % After this command, the following variables are present and
  % will be useful in the finite volume implementation.
  %
  % Nt: Number of triangles (i.e. elements) in mesh
  %
  % Nv: Number of nodes (i.e. vertices) in mesh
  %
  % Ne: Number of internal edges (between two triangles)
  %
  % Nbe: Number of edges which lie on a boundary of the computational
  %      domain.
  % 
  % tri2nod(3,Nt):  list of the 3 node numbers which form the current
  %                 triangle.  Thus, tri2nod(1,i) is the 1st node of
  %                 the i'th triangle, tri2nod(2,i) is the 2nd node
  %                 of the i'th triangle, etc.
  %
  % xy(2,Nv): list of the x and y locations of each node.  Thus,
  %           xy(1,i) is the x-location of the i'th node, xy(2,i)
  %           is the y-location of the i'th node, etc.
  %
  % edge2tri(1, i) = index of first  node of edge i
  % edge2tri(2, i) = index of second node of edge i
  % edge2tri(3, i) = index of triangle to left of edge i
  % edge2tri(4, i) = index of triangle to right of edge i
  %
  % The edge2tri nodes are oriented so that the triangle in 
  % edge2tri(3,i) is to the left when moving from node 1 to node 2
  %
  % bedge2tri(1, i) = index of first node for bedge i
  % bedge2tri(2, i) = index of second node for bedge i
  % bedge2tri(3, i) = index of the triangle for bedge i
  % bedge2tri(4, i) = boundary ID for bedge i 
  %                   (0 for cylinder, 1 for outer boundary)
  % 
  % The bedge2tri nodes are oriented so that the triangle is to 
  % the left when moving from node 1 to node 2

  % Enter freestream and initial Mach number
  Minf =  2.0;
  Minit = 0.0;
  
  % Specific heat ratio
  gamma = 1.4;
  
  % Initialize solution
  Uinit = [1; Minit; 0; 1/(gamma-1)/gamma + 0.5*Minit^2]; 
  U = ones(4,Nt);
  U(1,:) = Uinit(1)*U(1,:);
  U(2,:) = Uinit(2)*U(2,:);
  U(3,:) = Uinit(3)*U(3,:);
  U(4,:) = Uinit(4)*U(4,:);

  n = 0;

  % Clear figures
  figure(1); clf;
  figure(2); clf;

end

% Set tolerances on iteration and residual tolerance
if (~exist('ntol')),
  ntol = 100;
end
if (~exist('Rtol')),
  Rtol = -5;
end
Rmax = 0;

% Set CFL
if (~exist('CFL')),
  CFL = 0.9;
end

% Calculate interior edge unit normals and lengths

% INSERT CODE HERE.  

% Calculate boundary edge unit normals and lengths

% INSERT CODE HERE.  NOTE: The CalcForces.m routine requires
% that the boundary edge normals be stored in the array
% bnormal(i,j) where i = 1:2, and j = 1:Nbe.  The normals must
% point into the domain from the boundary.  Also, the boundary
% edge lengths must be stored in blength(j) where j = 1:Nbe.  

% Calculate cell areas

% INSERT CODE HERE NOTE: to use the cyl_adaptmesh routine,
% the cell areas must be stored in an array called A
% which is of length 1 x Nt (i.e. a row vector of length Nt).  
% This array has been initialized to zero in the above line.

% Calculate freestream states for use in setting boundary conditions
Uinf = [1; Minf; 0; 1/(gamma-1)/gamma + 0.5*Minf^2]; 

% Loop until one of the termination criterion are met
while ((n < ntol) & (Rmax > Rtol)),
  
  % Zero residual and dt
  R  = zeros(4,Nt);
  dt = zeros(1,Nt);
  
  % Calculate flux and smax at each interior edge 
  % sending flux to residuals and smax to dt
  for i = 1:Ne,

    iL = edge2tri(3,i); % Index of left triangle
    iR = edge2tri(4,i); % Index of right triangle
    
    % INSERT CODE HERE
    % Note: use the eulerflux function to calculate the flux
    % and smax
    
  end
  
  % Calculate flux and smax at each boundary edge 
  % sending flux to residuals and smax to dt
  for i = 1:Nbe,
    iL = bedge2tri(3,i); % Interior triangle
    if (bedge2tri(4,i) == 1), % Outer boundary
      % INSERT CODE HERE (call eulerflux)
    else, % Cylinder boundary
      % INSERT CODE HERE (call wallflux)
    end

    % INSERT CODE HERE TO UPDATE RESIDUAL and DT for
    % interior triangle

  end
  
  % Local timestep calculation (this line uses the
  % CFL definition and the accumulated values of smax
  % to calculate for each cell dt/A.
  
  % INSERT CODE HERE
  
  % Calculate maximum residual
  Rmax = log10(max(max(abs(R))));
  
  % Forward Euler step
  % INSERT CODE HERE

  % Increment iteration
  n = n + 1;

  % Calculate some useful quantities (Mach number)
  q = sqrt(U(2,:).^2 + U(3,:).^2)./U(1,:);
  p = (gamma-1)*(U(4,:) - 0.5*U(1,:).*q.^2);
  c = sqrt(gamma*p./U(1,:));
  M = q./c;
  
  % Calculate lift and drag coefficients
  CalcForces;
  
  % Plot current solution
  figure(1);
  clf;
  patch('Vertices',xy','Faces',tri2nod','CData',M','FaceColor', ...
        'flat','EdgeColor','none');
  axis equal;
  title(sprintf('Minf = %5.2f, CD = %4.2f, Log10 Max Res = %7.3f\n', ...
                Minf,CD,Rmax));
  xlabel('x'); ylabel('y');
  colorbar;
  
  % Plot max residual and drag coefficient
  figure(2);
  subplot(211);
  hold on; 
  plot(n,Rmax,'*'); 
  hold off;
  xlabel('Iteration');
  ylabel('Log10 Rmax');

  subplot(212);
  hold on;
  plot(n,CD,'*');
  hold off;
  xlabel('Iteration');
  ylabel('C_D');

  drawnow;

end

