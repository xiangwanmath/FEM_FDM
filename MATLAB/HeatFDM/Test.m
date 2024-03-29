tSteps = 10; % Number of iterations -> will replace this later
xSteps = 50; % Number of INTERIOR (BC non-inclusive) x coords
ySteps = 50; % Number of INTERIOR (BC non-inclusive) y coords
k = .00001; % Diffusion Constant 
tStep = .01; % dt
BCtype = ['d','d','d','d']; % Type of Boundary Conditions, d = dichirlet, n = neumann, r = robin
BoundCoord = [pi/2,(3*pi)/2; pi/2, (3*pi)/2]; % Top: X0->XL Bottom: Y0->YL

dif = Diffusion2D(tStep, tSteps, xSteps, ySteps, k, BoundCoord, BCtype);

d = dif.runSim();

%surf(d(:,:, 9));

X = dif.getXGrid();
Y = dif.getYGrid();

plt(X, Y, d, tSteps);


function plt(X, Y, TP, tIter)

for n=1:tIter

    surf(X, Y, TP(:, :, n));

    Mov(n)=getframe;

end

movie(Mov,1,1);

end


