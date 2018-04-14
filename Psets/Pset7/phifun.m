function phi = phifun(xj,dx)
%basis function
%   xj - xj
%   dx - change in nodal x

n = length(xj);
phi = cell(1,n);
points = n*49+1;
x=linspace(0,1,points);

for i = 1:n
    phi{i} = zeros(1,length(x));
    if i == 1
        ind = (x <= xj(i+1));
        phi{i} = (xj(i+1)-x)/dx;
        phi{i} = phi{i}.*ind;
    elseif i == n
        ind = (x >= xj(i-1));
        phi{i} = (x-xj(i-1))/dx;
        phi{i} = phi{i}.*ind;
    else
        for k = 1:length(x)
            if x(k)>=xj(i-1) && x(k)<=xj(i)
                phi{i}(k) = (x(k) - xj(i-1))/dx;
            elseif x(k)>xj(i) && x(k)<=xj(i+1)
                phi{i}(k) = (xj(i+1)-x(k))/dx;
            else
                phi{i}(k) = 0;
            end
        end
    end
end


end

