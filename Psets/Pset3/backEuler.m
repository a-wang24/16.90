function [x,y] = backEuler(xinit,yinit,xfinal,deltax)

    %f - this is your y prime
    %xinit - initial X
    %yinit - initial Y
    %xfinal - final X
    %deltat - step size

    n = (xfinal-xinit)/deltax; %Calculate steps

    %Inititialize arrays...
    %The first elements take xinit and yinit corespondigly, the rest fill with 0s.
    x = [xinit zeros(1,n)];
    y = [yinit zeros(1,n)];

    %Numeric routine
    for i = 1:n
        x(i+1) = x(i)+deltax;
        ynew = y(i)+deltax*(f(x(i),y(i)));
        y(i+1) = y(i)+deltax*f(x(i+1),ynew);
    end
end