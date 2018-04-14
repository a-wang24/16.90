%Alan Wang
%Pset 4 Script

%setup
%sample Nx and Nt
NxChoices = [10 20 40 80];
NtChoices = [50 100 200 400 800 1600 3200];
error = zeros(length(NtChoices),length(NxChoices));

%parameters
L = .2; %[cm]
rho = 2700; %[g/cm^3]
c = 900; %[J/gC]
k = 167; %[w/cmC]

for j=1:length(NtChoices)
    Nt = NtChoices(j);
    for k=1:length(NxChoices)
        
        %discretize x
        Nx = NxChoices(k);
        x = linspace(0,L,Nx+1);
        deltax = L/Nx;
        
        %Initial Condition
        Tx0 = 20 + 100*sin(pi*x/L);
        
        t = 0;
        tcool = log(10)*((rho*c*L^2)/(k*pi^2));
        tfinal = tcool;
        deltat = tfinal/Nt;
        T = Tx0;
        Txt = 20 + 100*exp(((-k*pi^2)/(rho*c*L^2))*tcool)*sin(pi*x/L);
        
        while t<tfinal
            T(1) = 20;
            T(end) = 20;
            T(2:end-1) = T(2:end-1) + ((deltat*k)/(rho*c*deltax^2))*(T(3:end)-2*T(2:end-1)+T(1:end-2));
            t = t + deltat;
        end
        
        Terr = abs(T-Txt);
        tempErr = 0;
        
        for i = 1:(Nx+1)
            tempErr = tempErr + Terr(i);
        end
        
        tempErr = (tempErr)/(Nx+1);
        error(j,k) = tempErr;
    end
end

for i = 1:length(NtChoices)
    for j = 1:length(NxChoices)
        disp(['Nt=',num2str(NtChoices(i)),' Nx=',num2str(NxChoices(j)),' error= ',num2str(error(i,j))])
    end
end
