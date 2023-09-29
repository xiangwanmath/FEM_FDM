gridSize = 50;

steps = [gridSize, gridSize];
lengths = [8*pi, 8*pi];

test = Laplace2DFD(steps, lengths, lengths(1)/steps(1));

M = test.runSimNeu();
surf(M);
% Look up sparse matrices for effiency 

