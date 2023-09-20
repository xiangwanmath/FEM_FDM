gridSize = 50;

steps = [gridSize, gridSize];
lengths = [4*pi, 4*pi];

test = Laplace2DFD(steps, lengths, lengths(1)/steps(1));

M = test.runSim();

surf(M);

% Look up sparse matrices for effiency 
% Log(err) = Log(h)
% Find regression
