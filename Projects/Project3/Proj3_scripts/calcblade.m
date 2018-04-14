function [T, Tgrad] = calcblade(params, chord, tri2nod, xy, bedge)

% Set blade material properties
kblade = params(1);

% Set gas temperature and wall heat transfer coefficients at
% boundaries of the blade.  Note: Tcool(i) and hcool(i) are the
% values of Tcool and hcool for the ith boundary which are numbered
% as follows:  
%
%   1 = 1st internal cooling passage (from leading edge)
%   2 = 2nd internal cooling passage 
%   3 = 3rd internal cooling passage (including trailing edge slot)

Tgas   = params(2);
hgasLE = params(3);
hgasTE = params(4);

Tcool(1:3) = params(5); % Values assumed to be the same for all passages
hcool(1:3) = params(6); % Values assumed to be the same for all passages

% Find problem size
[Ntmp, Nv]  = size(xy);
[Ntmp, Nt]  = size(tri2nod);
[Ntmp, Nbc] = size(bedge);

% Zero stiffness matrix

K = sparse(Nv, Nv);
b = zeros(Nv, 1);

% Loop over elements and calculate residual and stiffness matrix

for ii = 1:Nt,
  
  kn(1) = tri2nod(1,ii);
  kn(2) = tri2nod(2,ii);
  kn(3) = tri2nod(3,ii);
    
  xe(1) = xy(1,kn(1));
  xe(2) = xy(1,kn(2));
  xe(3) = xy(1,kn(3));

  ye(1) = xy(2,kn(1));
  ye(2) = xy(2,kn(2));
  ye(3) = xy(2,kn(3));

  % Calculate all of the necessary shape function derivatives, the
  % Jacobian of the element, etc.
  
  % Derivatives of node 1's interpolant 
  dNdxi(1,1) = -1.0; % with respect to xi1
  dNdxi(1,2) = -1.0; % with respect to xi2
  
  % Derivatives of node 2's interpolant
  dNdxi(2,1) =  1.0; % with respect to xi1
  dNdxi(2,2) =  0.0; % with respect to xi2

  % Derivatives of node 3's interpolant
  dNdxi(3,1) =  0.0; % with respect to xi1
  dNdxi(3,2) =  1.0; % with respect to xi2
  
  % Sum these to find dxdxi (note: these are constant within an element)
  dxdxi = zeros(2,2);
  for nn = 1:3,
    dxdxi(1,:) = dxdxi(1,:) + xe(nn)*dNdxi(nn,:);
    dxdxi(2,:) = dxdxi(2,:) + ye(nn)*dNdxi(nn,:);
  end
  
  % Calculate determinant for area weighting
  J = abs(dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1));
  A = 0.5*J; % Area is half of the Jacobian
  
  % Invert dxdxi to find dxidx using inversion rule for a 2x2 matrix
  dxidx = [ dxdxi(2,2)/J, -dxdxi(1,2)/J; ...
	   -dxdxi(2,1)/J,  dxdxi(1,1)/J];
  
  % Calculate dNdx 
  dNdx = dNdxi*dxidx;

  % Add contributions to stiffness matrix for node 1 weighted residual
  K(kn(1), kn(1)) = K(kn(1), kn(1)) - kblade*(dNdx(1,1)*dNdx(1,1) + dNdx(1,2)*dNdx(1,2))*A;
  K(kn(1), kn(2)) = K(kn(1), kn(2)) - kblade*(dNdx(1,1)*dNdx(2,1) + dNdx(1,2)*dNdx(2,2))*A;
  K(kn(1), kn(3)) = K(kn(1), kn(3)) - kblade*(dNdx(1,1)*dNdx(3,1) + dNdx(1,2)*dNdx(3,2))*A;
  
  % Add contributions to stiffness matrix for node 2 weighted residual
  K(kn(2), kn(1)) = K(kn(2), kn(1)) - kblade*(dNdx(2,1)*dNdx(1,1) + dNdx(2,2)*dNdx(1,2))*A;
  K(kn(2), kn(2)) = K(kn(2), kn(2)) - kblade*(dNdx(2,1)*dNdx(2,1) + dNdx(2,2)*dNdx(2,2))*A;
  K(kn(2), kn(3)) = K(kn(2), kn(3)) - kblade*(dNdx(2,1)*dNdx(3,1) + dNdx(2,2)*dNdx(3,2))*A;
  
  % Add contributions to stiffness matrix for node 3 weighted residual
  K(kn(3), kn(1)) = K(kn(3), kn(1)) - kblade*(dNdx(3,1)*dNdx(1,1) + dNdx(3,2)*dNdx(1,2))*A;
  K(kn(3), kn(2)) = K(kn(3), kn(2)) - kblade*(dNdx(3,1)*dNdx(2,1) + dNdx(3,2)*dNdx(2,2))*A;
  K(kn(3), kn(3)) = K(kn(3), kn(3)) - kblade*(dNdx(3,1)*dNdx(3,1) + dNdx(3,2)*dNdx(3,2))*A;

end

% Loop over boundary edges and account for bc's
% Note: the bc's are all convective heat transfer coefficient bc's
% so the are of 'Robin' form.  This requires modification of the
% stiffness matrix as well as impacting the right-hand side, b.
%
for ii = 1:Nbc,
  
  % Get node numbers on edge
  kn(1) = bedge(1,ii);
  kn(2) = bedge(2,ii);
  
  % Get node coordinates
  xe(1) = xy(1,kn(1));
  xe(2) = xy(1,kn(2));
  
  ye(1) = xy(2,kn(1));
  ye(2) = xy(2,kn(2));
  
  % Calculate edge length
  ds = sqrt((xe(1)-xe(2))^2 + (ye(1)-ye(2))^2);
  
  % Determine the boundary number
  bnum = bedge(3,ii);
  if (bnum == 0), % Blade surface
    
    [Text, hext] = Thgas( 0.5*(xe(1)+xe(2)), 0.5*(ye(1)+ye(2)), Tgas, ...
                          hgasLE, hgasTE, chord );
    
  else, % cooling passage
    Text = Tcool(bnum);
    hext = hcool(bnum);
    
  end
  
  % Based on boundary number, set heat transfer bc
  K(kn(1), kn(1)) = K(kn(1), kn(1)) - hext*ds*(1/3);
  K(kn(1), kn(2)) = K(kn(1), kn(2)) - hext*ds*(1/6);
  b(kn(1))        = b(kn(1))        + hext*ds*0.5*Text;
  
  K(kn(2), kn(1)) = K(kn(2), kn(1)) - hext*ds*(1/6);
  K(kn(2), kn(2)) = K(kn(2), kn(2)) - hext*ds*(1/3);
  b(kn(2))        = b(kn(2))        + hext*ds*0.5*Text;
  
end

% Solve for temperature
T = -K\b;

% Post-process to find temperature gradients in each element
Tgrad = zeros(Nt, 1);
for ii = 1:Nt,
  
  kn(1) = tri2nod(1,ii);
  kn(2) = tri2nod(2,ii);
  kn(3) = tri2nod(3,ii);
    
  xe(1) = xy(1,kn(1));
  xe(2) = xy(1,kn(2));
  xe(3) = xy(1,kn(3));

  ye(1) = xy(2,kn(1));
  ye(2) = xy(2,kn(2));
  ye(3) = xy(2,kn(3));

  % Calculate all of the necessary shape function derivatives, the
  % Jacobian of the element, etc.
  
  % Derivatives of node 1's interpolant 
  dNdxi(1,1) = -1.0; % with respect to xi1
  dNdxi(1,2) = -1.0; % with respect to xi2
  
  % Derivatives of node 2's interpolant
  dNdxi(2,1) =  1.0; % with respect to xi1
  dNdxi(2,2) =  0.0; % with respect to xi2

  % Derivatives of node 3's interpolant
  dNdxi(3,1) =  0.0; % with respect to xi1
  dNdxi(3,2) =  1.0; % with respect to xi2
  
  % Sum these to find dxdxi (note: these are constant within an element)
  dxdxi = zeros(2,2);
  for nn = 1:3,
    dxdxi(1,:) = dxdxi(1,:) + xe(nn)*dNdxi(nn,:);
    dxdxi(2,:) = dxdxi(2,:) + ye(nn)*dNdxi(nn,:);
  end
  
  % Calculate determinant for area weighting
  J = abs(dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1));
  A = 0.5*J; % Area is half of the Jacobian
  
  % Invert dxdxi to find dxidx using inversion rule for a 2x2 matrix
  dxidx = [ dxdxi(2,2)/J, -dxdxi(1,2)/J; ...
	   -dxdxi(2,1)/J,  dxdxi(1,1)/J];
  
  % Calculate dNdx 
  dNdx = dNdxi*dxidx;

  % Calculate x and y derivatives
  Tx = T(kn(1))*dNdx(1,1) + T(kn(2))*dNdx(2,1) + T(kn(3))*dNdx(3,1);
  Ty = T(kn(1))*dNdx(1,2) + T(kn(2))*dNdx(2,2) + T(kn(3))*dNdx(3,2);

  % Calculate gradient magnitude
  Tgrad(ii) = sqrt(Tx*Tx + Ty*Ty);
  
end

