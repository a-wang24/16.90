function fu = fufun(Un)
%ffun represents the f(u) of the problem
%   Un is a state vector
global M_hh M_ha M_aa M_ah D_h D_a K_h K_a K_nl

a = Un(1);
h = Un(2);
p = Un(3);
v = Un(4);
L = @Lfun;
M = @Mfun;
fu = [p; 
    v; 
    (M_ah*D_h*v + M_ah*K_h*h + M_ah*L(Un) - M_hh*D_a*p - K_a*M_hh*(1+K_nl*(h^2))*a - M(Un)*M_hh)/(M_aa*M_hh - M_ah*M_ha);
    (M_aa*D_h*v + M_aa*K_h*h + M_aa*L(Un) - M_ha*D_a*p - K_a*M_ha*(1+K_nl*(h^2))*a - M(Un)*M_ha)/(M_ah*M_ha - M_aa*M_hh)
    ];

end

