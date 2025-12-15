% ----- Parameters -----

% Number of grid points on the inflow boundary (the grid points in the 
% rest of the domain is scaled accordingly)
m = 15;

% Reynolds number
Re = 80;

% End time 
Tend = 6;

% Length of outflow pipe %från början 2
Lout = 5;

% Show animation and streamlines of the solution? True or false. 
animate = true;

% ----- Main function -----
incomp(m,Re,Tend,Lout,animate)