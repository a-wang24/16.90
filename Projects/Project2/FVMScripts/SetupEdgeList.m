function [edge2tri, bedge2tri] = SetupEdgeList(tri2nod, bedge)

% This routine determines the edge data structures from a
% triangular mesh connectivity pattern (stored in tri2nod)
% and a list of boundary edge nodes (stored in bedge).  The
% resulting data structures are edge2tri and bedge2tri which store
% the following info:
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
% 
% The bedge2tri nodes are oriented so that the triangle is to the left
% when moving from node 1 to node 2
%
% Note the edges in bedge2tri are stored in the same order as they
% appear in the input array bedge2tri
%

% Get sizes
[Ntemp, Nt] = size(tri2nod);
[Ntemp, Nbe] = size(bedge);

% Determine number of internal edges
Ne = (3*Nt - Nbe)/2;

% Determine maximum node number
Nv = max(max(tri2nod));

% Set up arrays to hold edge lists
edge2tri  = zeros(4,Ne);
bedge2tri = zeros(4,Nbe);

% Copy over contents of bedge to bedge2tri
bedge2tri(1:2,:) = bedge(1:2,:);
bedge2tri(4,:)   = bedge(3,:);

% Set up work space 
nod2edge = sparse(Nv,Nv);

% Set interior triangle index in bedge2tri list
for i = 1:Nbe,

  % Get node numbers
  k1 = bedge2tri(1,i);
  k2 = bedge2tri(2,i);
  
  % Find min and max node numbers
  kmin = min(k1, k2);
  kmax = max(k1, k2);
  
  % Check if edge does not exist in nod2edge list
  ktri = nod2edge(kmin,kmax);
  if (ktri ~= 0), % Edge exists: Error
    fprintf('Multiple boundary edge. ERROR: ktri ~= 0.  ktri = %i\n',ktri);
    return;
  else, % Edge does not exist, insert it
    nod2edge(kmin,kmax) = -i; % Add edge to list noting triangle
  end
end

% Loop over triangles and find all edges
ie = 0;
for i = 1:Nt,
  for j = 1:3,
    
    % Get node numbers
    k1 = tri2nod(         j, i);
    k2 = tri2nod(mod(j,3)+1, i);
    
    % Find min and max node numbers
    kmin = min(k1, k2);
    kmax = max(k1, k2);
    
    % Check if edge does not exist in nod2edge list
    ktri = nod2edge(kmin,kmax);
    if (ktri == 0), % Edge does not exist
      nod2edge(kmin,kmax) = i; % Add edge to list noting triangle
    elseif (ktri > 0), % Edge exists, insert into edge2tri list
      ie = ie + 1;
      edge2tri(1,ie) = kmin;
      edge2tri(2,ie) = kmax;
      if (kmin == k1),
        edge2tri(3,ie) = i;
        edge2tri(4,ie) = ktri;
      else
        edge2tri(3,ie) = ktri;
        edge2tri(4,ie) = i;
      end
      nod2edge(kmin,kmax) = 0; 
    else, % Boundary edge exists, insert into bedge2tri list
      bedge2tri(1,-ktri) = k1;
      bedge2tri(2,-ktri) = k2;
      bedge2tri(3,-ktri) = i;
      nod2edge(kmin,kmax) = 0;
    end
    
  end
end

% Check that no non-zero entries exist in nod2edge
nz = nnz(nod2edge);
if (nz > 0),
  fprintf('ERROR: non-zero entries remain in nod2edge\n');
end

