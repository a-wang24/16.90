function R = Resfun(wm,vn,vn1,deltax,dfdx)
%Resfun calculates Residual for 2nd Order Backward Differentiation
%wm     - current guess
%vn     - latest solution
%vn1    - previous solution
%deltax - time step
%dfdx   - function handle for df/dx

R = wm - (4/3)*vn + (1/3)*vn1 - 2/3*deltax*dfdx(wm);

end

