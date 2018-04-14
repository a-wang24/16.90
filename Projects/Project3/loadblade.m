function [tri2nod, xy, bedge] = loadblade(fname, chord)

% Load in the grid file
load(fname);

% NOTE:  after loading a gridfile using the load(fname) command,
%        three important grid variables and data arrays exist.  These are:
%
% tri2nod(3,Nt):  list of the 3 node numbers which form the current
%                 triangle.  Thus, tri2nod(1,i) is the 1st node of
%                 the i'th triangle, tri2nod(2,i) is the 2nd node
%                 of the i'th triangle, etc.
%
% xy(2,Nv): list of the x and y locations of each node.  Thus,
%           xy(1,i) is the x-location of the i'th node, xy(2,i)
%           is the y-location of the i'th node, etc.
%
% bedge(3,Nbc): For each boundary edge, bedge(1,i) and bedge(2,i) 
%               are the node numbers for the nodes at the end
%               points of the i'th boundary edge.  bedge(3,i) is an
%               integer which identifies which boundary the edge is
%               on. In this solver, the third value has the
%               following meaning:
%
%               bedge(3,i) = 0: edge is on the airfoil
%               bedge(3,i) = 1: edge is on the first cooling passage
%               bedge(3,i) = 2: edge is on the second cooling passage
%               bedge(3,i) = 3: edge is on the third cooling passage
% 

% Re-scale xy location because they currently have unit chord.
xy = chord*xy;

