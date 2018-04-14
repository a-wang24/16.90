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

%% Task 1 - Estimating Probability of Failure
% Run a Monte Carlo Simulation
% See MonteCarloWhile.m for while loop version without pre calculated N

% Create vectors to hold input parameters
N = 13000;
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

for i = 1:N
    % Set input parameter values
    % Gives them a sampled value based on their triangular distribution
    % Each iteration will use different input value
    kblade(i) =    trd(28.5,30,31.5);
    Tgas(i)   =  trd(1400,1500,1600);
    hgasLE(i) = trd(15500,16000,18000);
    hgasTE(i) =  trd(3000,4000,5000);
    Tcool(i)  =   trd(550,600,650);
    hcool(i)  =  trd(1400,1500,1600);
    
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
end

%% Post Process Simulation Data for probability estimates
% Calculate probability of failure estimate
probFail = failcount/N;

% Find estimated probability for Tmax > Tlimit and gradTmax > gradTl
Tcutoff = Tmax > Tl;
Tcount = nnz(Tcutoff);
Tprob = Tcount/N;
gradTcutoff = Tgradmax > gradTl;
gradTcount = nnz(gradTcutoff);
gradTprob = gradTcount/N;

% Calculate 99% confidence intervals for probability of Tmax > Tlimit
% and probability of gradTmax > gradTlimit
% Here we use that estimate for probability is normal with mean actual
% probability of event occurring and standard error sqrt(P(E)*(1-P(E))/N)
TmgC1 = Tprob - 3*sqrt((Tprob)*(1-Tprob)/N);
TmgC2 = Tprob + 3*sqrt((Tprob)*(1-Tprob)/N);
gradTmgC1 = gradTprob - 3*sqrt((gradTprob)*(1-gradTprob)/N);
gradTmgC2 = gradTprob + 3*sqrt((gradTprob)*(1-gradTprob)/N);

%% Plot histograms and CDF plots for Tmax, gradTmax and Damage

figure(1)
histogram(Tmax,50, 'Normalization','pdf'); %histogram with 50 bins of Tmax
xlabel('Tmax [K]')

figure(2)
histogram(Tgradmax,50, 'Normalization','pdf');
xlabel('Tgradmax [K/m]')

figure(3)
histogram(D,50, 'Normalization','pdf');
xlabel('Damage')

figure(4)
subplot(3,1,1)
cdfplot(Tmax)
title('CDF of Tmax')
subplot(3,1,2)
cdfplot(Tgradmax)
title('CDF of TgradMax')
subplot(3,1,3)
cdfplot(D)
title('CDF of Damage')

%% Task 2 - Screening for Important Factors
% Calculate Spearman Correlation Coefficients

% sort damage and input parameters from low to high
% ID is vector of ranks, for example ID(1) = rank(D(1))
[sortD,ID] = sort(D);
[sortKblade,IK] = sort(kblade);
[sortTgas, ITg] = sort(Tgas);
[sorthgasLE, IhgL] = sort(hgasLE);
[sorthgasTE, IhgT] = sort(hgasTE);
[sortTcool, ITc] = sort(Tcool);
[sorthcool, Ihc] = sort(hcool);

% compute correlation coefficients between damage and each input
% correlation coefficient between damage and kblade
dDK = ID - IK; %calculate rank distance for each pair
dDKsquared = dDK.^2; %square each element value
sumDK = sum(dDKsquared); %sum them all up
rhoDK = 1 - (6*sumDK)/(N*((N^2) - 1)); %use value in correlation coeff formula

% correlation coefficient between damage and Tgas
dDTg = ID - ITg;
dDTgsquared = dDTg.^2;
sumDTg = sum(dDTgsquared);
rhoDTg = 1 - (6*sumDTg)/(N*((N^2) - 1));

% correlation coefficient between damage and hgasLE
dDhgL = ID - IhgL;
dDhgLsquared = dDhgL.^2;
sumDhgL = sum(dDhgLsquared);
rhoDhgL = 1 - (6*sumDhgL)/(N*((N^2) - 1));

% correlation coefficient between damage and hgasTE
dDhgT = ID - IhgT;
dDhgTsquared = dDhgT.^2;
sumDhgT = sum(dDhgTsquared);
rhoDhgT = 1 - (6*sumDhgT)/(N*((N^2) - 1));

% correlation coefficient between damage and Tcool
dDTc = ID - ITc;
dDTcsquared = dDTc.^2;
sumDTc = sum(dDTcsquared);
rhoDTc = 1 - (6*sumDTc)/(N*((N^2) - 1));

% correlation coefficient between damage and hcool
dDhc = ID - Ihc;
dDhcsquared = dDhc.^2;
sumDhc = sum(dDhcsquared);
rhoDhc = 1 - (6*sumDhc)/(N*((N^2) - 1));

%% Task 3 - Impact of Important Factors on Probability of Failure
% reduce input variability and run MC simulation

% Create vectors to hold input parameters
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

for i = 1:N
    % Set input parameter values
    % Gives them a sampled value based on their triangular distribution
    % Each iteration will use different input value
    % Change the min, mpp, max values below accordingly
    kblade(i) =    trd(28.5,30,31.5);
    Tgas(i)   =  trd(1400,1500,1600);
    hgasLE(i) = trd(15500,16000,18000);
    hgasTE(i) =  trd(3000,4000,5000);
    Tcool(i)  =   trd(550,600,650);
    hcool(i)  =  trd(1450,1500,1550);
    
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
    % accounts for +- .01 error with 99% confidence
    % breaks when 3 standard deviations away is less than .01
    if i > 100 && standErr <= 1/300
        break;
    end
end

probFail = failcount/N;