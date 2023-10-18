function cm = redwhiteblue
% A colormap going from red to blue via white in 256 steps.
%

% Number of points
n = 128;

% Generate red to white
R1 = linspace(1, 1, n)';
G1 = linspace(0, 1, n)';
B1 = linspace(0, 1, n)';

% Generate white to blue
R2 = linspace(1, 0, n)';
G2 = linspace(1, 0, n)';
B2 = linspace(1, 1, n)';

% Concatenate RGB vectors
R = [R1; R2];
G = [G1; G2];
B = [B1; B2];

% Create colormap
cm = [R, G, B];
