gridSize = 70;
steps = [gridSize, gridSize];
BCoords = [pi/2, (3*pi)/2; pi/2, (3*pi)/2];
BC = "r";
BCtype = [BC,BC,BC,BC];

test = Laplace2DFD(steps, BCoords, BCtype);

A = test.makeA();

M = test.runSim();
surf(M);

