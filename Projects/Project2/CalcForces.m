% Calculate force coefficients on cylinder

CD = 0;
CL = 0;
for i = 1:Nbe,
  iL = bedge2tri(3,i);
  if (bedge2tri(4,i) == 0), % Cylinder boundary
    [H, smax] = wallflux(U(:,iL),bnormal(:,i),gamma);
    CD = CD - H(2)*blength(i);
    CL = CL - H(3)*blength(i);
  end
end
CD = CD/(Minf^2);
CL = CL/(Minf^2);
