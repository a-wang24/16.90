function [x] = trirnd(xmin,xmpp,xmax)

% Generate a random number from a triangular distribution
% in which:
%
% xmin = minimum value of x
% xmpp = most-probable value of x
% xmax = maximum value of x
%

% Find percentile from uniform distribution
u = rand;

% Calculate percentile at xmpp
Fmpp = (xmpp-xmin)/(xmax-xmin);

% Check is u < Fmpp, and find x = F(u)
if (u < Fmpp), 
  x = sqrt(u*(xmax-xmin)*(xmpp-xmin)) + xmin;
else
  x = xmax - sqrt((1-u)*(xmax-xmin)*(xmax-xmpp));
end


