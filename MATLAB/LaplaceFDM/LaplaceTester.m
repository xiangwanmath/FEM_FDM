

gridSize = 50;
steps = [gridSize, gridSize];
BCoords = [pi/2, (3*pi)/2; pi/2, (3*pi)/2];
BC = "d";
BCtype = [BC,BC,BC,BC];

test = Laplace2DFD(steps, BCoords, BCtype);

A = test.makeA();
%disp(A);
M = test.runSim();
surf(M);



%{
radius = 1;
steps = [4, 4]; % 1 is r 2 is t

test = LaplacePOLAR(radius,steps,"d");
U = test.runSim();
disp(U);
%}