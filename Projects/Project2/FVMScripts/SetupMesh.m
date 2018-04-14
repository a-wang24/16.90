Nv  = max(size(cyl_p));
Nbc = max(size(cyl_e));
Nt  = max(size(cyl_t));

% Put it into the expected form
xy = cyl_p;

tri2nod = cyl_t(1:3,:);

bedge = [cyl_e(1:2,:); zeros(1,Nbc)];

for ii = 1:Nbc,
  rb = norm(xy(:,bedge(1,ii))); % Calculate radius of a boundary point
  if (rb > 1.5), % Assume outer boundary is at least 0.5 away from cylinder
    bedge(3,ii) = 1; % Outer boundary
  end
end

% Form edge list for use in finite volume method
[edge2tri, bedge2tri] = SetupEdgeList(tri2nod, bedge);
Ne  = length(edge2tri);
Nbe = length(bedge2tri);
