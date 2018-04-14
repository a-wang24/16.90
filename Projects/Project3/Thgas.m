function [Tgas, hgas] = Thgas(x, y, Tgas, hgasLE, hgasTE, chord)

% This function calculate the gas path temperature
% and heat transfer coefficient given the x and y location.

xc = x/chord;

% Length of high heat transfer region at leading-edge
lhot = 0.05;

hgas = hgasTE + (hgasLE-hgasTE)*exp(-4*(xc/lhot)^2);

