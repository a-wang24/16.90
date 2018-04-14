%Alan Wang
%16.90 Project 3

%%Setup
% Setup the varability of input parameters with triangular distributions
% and other parameters
trd = @trirnd;
Tl = 1500; %[K]
gradTl = 80000; %[K/m]
Tdm = 1430; %[K]
gradTdm = 70000; %[K/m]

% Load blade geometry and mesh
chord = 0.04;
[tri2nod, xy, bedge] = loadblade('hpblade_coarse', chord);

N = 20000;
Tmax = zeros(1,N);
Tgradmax = zeros(1,N);
D = zeros(1,N);
kblade = zeros(1,N);
Tgas = zeros(1,N);
hgasLE = zeros(1,N);
hgasTE = zeros(1,N);
Tcool = zeros(1,N);
hcool = zeros(1,N);
failcount = 0;
i = 1;
standErr = 0;

while i < 100 || standErr >= 1/300
    % Set input parameter values
    % Gives them a sampled value based on their triangular distribution
    % Each iteration will use different input value
    kblade(i) =    trd(28.5,30,31.5);
    Tgas(i)   =  trd(1400,1500,1600);
    hgasLE(i) = trd(15500,16000,18000);
    hgasTE(i) =  trd(3000,4000,5000);
    Tcool(i)  =   trd(550,600,650);
    hcool(i)  =  trd(1400,1500,1500);
    
    % Store parameters in a vector
    params = [kblade(i), Tgas(i), hgasLE(i), hgasTE(i), Tcool(i), hcool(i)];
    
    % Calculate temperature and gradient magnitude in blade (FEM analysis)
    [T, Tgrad] = calcblade( params, chord, tri2nod, xy, bedge);
    
    % Find max T and max gradient
    tempTmax = max(T);
    tempTgradmax = max(Tgrad);
    Tmax(i) = tempTmax;
    Tgradmax(i) = tempTgradmax;
    % Calculate Damage to blade
    tempD = ((tempTmax-Tdm)/(Tl - Tdm)) + ((tempTgradmax - gradTdm)/(gradTl - gradTdm));
    D(i) = tempD;
    % if Damage is above 1 then blade failed, increase event counter
    if tempD > 1
        failcount = failcount+1;
    end
    tempProb = failcount/i;
    standErr = sqrt((tempProb*(1-tempProb))/i);
    i = i+1;
    disp(i)
end

probFail = failcount/i;