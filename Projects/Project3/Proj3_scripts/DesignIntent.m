% Set input parameter values
kblade =    30;
Tgas   =  1500;
hgasLE = 16000;
hgasTE =  4000;
Tcool  =   600;
hcool  =  1500;

% Store parameters in a vector
params = [kblade, Tgas, hgasLE, hgasTE, Tcool, hcool];

% Load blade geometry and mesh
chord = 0.04;
[tri2nod, xy, bedge] = loadblade('hpblade_coarse', chord);

% Calculate temperature and gradient magnitude in blade (FEM analysis)
[T, Tgrad] = calcblade( params, chord, tri2nod, xy, bedge);

% Find max T and max gradient
Tmax = max(T);
Tgradmax = max(Tgrad);

fprintf('Tmax = %f, Tgradmax = %f\n',Tmax,Tgradmax);

