function dfdu = dfdufun(Un)
%Jacobian of f(u)
%   
global M_hh M_ha M_aa M_ah D_h D_a K_h K_a K_nl Q

a = Un(1);
h = Un(2);
cl = Q;
cm = -0.7*Q;
denom = (M_ah*M_ha - M_aa*M_hh);

dfdu = [0 0 1 0;
        0 0 0 1;
        (M_hh*(cm+K_a*(1+K_nl*h^2))-M_ah*cl)/denom (2*a*h*K_a*K_nl*M_hh-K_h*M_ah)/denom (D_a*M_hh)/denom (-D_h*M_ah)/denom;
        (-M_ha*(cm+K_a*(1+K_nl*h^2))+M_aa*cl)/denom (-2*a*h*K_a*K_nl*M_ha+K_h*M_aa)/denom (-D_a*M_ha)/denom (D_h*M_aa)/denom];
    
end

