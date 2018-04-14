% This Matlab script solves the two-dimensional Euler
% equations using a finite volume algorithm for the flow
% around a cylinder.
%

%set FVMrestart to 1 to iterate from current solution in memory

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
    ntol = 1000;
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

%Initiate variables for unit normals and lengths for interior edges
inormal = zeros(2,Ne);
ilengths = zeros(1,Ne);

%Calculate unit normals and lengths of interior edges and populate vectors
for i = 1:Ne
    %finds indices of nodes for edge i
    node1 = edge2tri(1,i);
    node2 = edge2tri(2,i);
    %finds x and y position of nodes for edge i
    x1 = xy(1,node1);
    y1 = xy(2,node1);
    x2 = xy(1,node2);
    y2 = xy(2,node2);
    %calculate lengths
    dist = sqrt(((x2-x1)^2)+((y2-y1)^2));
    ilengths(1,i) = dist;
    %calculate unit normal
    vect = [(x2-x1) (y2-y1)];
    normVect = [-vect(2) vect(1)];
    unitNormVect = normVect./dist;
    inormal(:,i) = unitNormVect;
end

% Calculate boundary edge unit normals and lengths

%Initiate variables for boundary lengths and normals
bnormal = zeros(2,Nbe);
blength = zeros(1,Nbe);

for i = 1:Nbe
    %finds indices of nodes for boundary edge i
    node1 = bedge2tri(1,i);
    node2 = bedge2tri(2,i);
    %finds x and y position of nodes for edge i
    x1 = xy(1,node1);
    y1 = xy(2,node1);
    x2 = xy(1,node2);
    y2 = xy(2,node2);
    %calculate lengths
    dist = sqrt(((x2-x1)^2)+((y2-y1)^2));
    blength(1,i) = dist;
    %calculate unit normal
    vect = [(x2-x1) (y2-y1)];
    normVect = [-vect(2) vect(1)];
    unitNormVect = normVect./dist;
    bnormal(:,i) = unitNormVect;
end

% Calculate cell areas

%initiate variable for area
A = zeros(1,Nt);

%populate the row vector A
for i = 1:Nt
    %find index of 3 nodes of triangle i
    node1 = tri2nod(1,i);
    node2 = tri2nod(2,i);
    node3 = tri2nod(3,i);
    %get x and y position of each node
    x1 = xy(1,node1);
    y1 = xy(2,node1);
    x2 = xy(1,node2);
    y2 = xy(2,node2);
    x3 = xy(1,node3);
    y3 = xy(2,node3);
    %calculate area of triangle i
    area = abs((x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2);
    A(1,i) = area;
end

% Calculate freestream states for use in setting boundary conditions
Uinf = [1; Minf; 0; 1/(gamma-1)/gamma + 0.5*Minf^2];

% Loop until one of the termination criterion are met
while ((n < ntol) && (Rmax > Rtol)),
    
    % Zero residual and dt
    R  = zeros(4,Nt);
    dt = zeros(1,Nt);
    
    %Zero smaxi*hi value for each cell
    sumSmaxH = zeros(1,Nt);
    
    % Calculate flux and smax at each interior edge
    % sending flux to residuals and smax to dt
    for i = 1:Ne,
        
        iL = edge2tri(3,i); % Index of left triangle
        iR = edge2tri(4,i); % Index of right triangle
        
        %establish state vectors for left and right triangle
        %find normal calculated above and use function eulerflux
        UL = U(:,iL);
        UR = U(:,iR);
        normal = inormal(:,i);
        [H,smax] = eulerflux(UL,UR,normal,gamma);
        %find h in ilength calculated above
        h = ilengths(1,i);
        %update Residual array with calculated flux
        R(:,iL) = R(:,iL) - H*h;
        R(:,iR) = R(:,iR) + H*h;
        %update max prop speed if less than previous max prop speed
        if smax > dt(1,iL)
            dt(1,iL) = smax;
        end
        %populate sum of smax*hi vector for later use when implementing
        %local timestepping
        sumSmaxH(1,iL) = sumSmaxH(1,iL) + smax*h;
        sumSmaxH(1,iR) = sumSmaxH(1,iR) + smax*h;
    end
    
    % Calculate flux and smax at each boundary edge
    % sending flux to residuals and smax to dt
    for i = 1:Nbe,
        iL = bedge2tri(3,i); % Interior triangle
        UL = U(:,iL); %state vector for left triangle
        normal = bnormal(:,i); %normal vector
        if (bedge2tri(4,i) == 1), % Outer boundary
            UR = Uinf; %right triangle state vector is freestream
            [H,smax] = eulerflux(UL,UR,normal,gamma);
        else % Cylinder boundary
            [H,smax] = wallflux(UL,normal,gamma);
        end
        
        %update residual and dt
        h = blength(1,i);
        R(:,iL) = R(:,iL) - H*h;
        if smax > dt(1,iL)
            dt(1,iL) = smax;
        end
        %again, find h and populate sum of smax*h for use when implementing
        %local timestepping
        sumSmaxH(1,iL) = sumSmaxH(1,iL) + smax*h;
    end
    
    % Local timestep calculation (this line uses the
    % CFL definition and the accumulated values of smax
    % to calculate for each cell dt/A.
    deltatRow = (2*CFL)./sumSmaxH;
    deltat = zeros(4,Nt);
    deltat(1,:) = deltatRow;
    deltat(2,:) = deltatRow;
    deltat(3,:) = deltatRow;
    deltat(4,:) = deltatRow;
    
    % Calculate maximum residual
    Rmax = log10(max(max(abs(R))));
    
    % Forward Euler step
    U = U - deltat.*R;
    
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

