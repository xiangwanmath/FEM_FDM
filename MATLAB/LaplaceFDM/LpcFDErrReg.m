length = [4*pi, 4*pi]; 
steps = [50, 50];
h = zeros(5, 1); % try to make h smaller
uEst = [];
uExact = [];


for i = 1:5
    h(i) = length(1)/steps(1);
    uEst = [uEst; uXYEst(steps(1), length(1), h(i))];
    steps = steps + 10;
end

steps = [50, 50];

for i = 1:5
    u = zeros(steps(1)^2, 1);
    for j = 1:steps(1)
        for k = 1:steps(1)
            u(j*k, 1) = uXY(j*h(i), k*h(i));
        end
    end
    uExact = [uExact; u];
    steps = steps + 10;
end

uErr = abs(uExact - uEst); 
uErr = log(uErr);

h = log(h);

steps = [50, 60, 70, 80, 90];

hV = [];

for i = 1:5
    hV = [hV; repmat(h(i), steps(i)^2, 1)];
end

slp = hV\uErr;

disp(slp);

function u = uXYEst(steps, length, h)
    est = Laplace2DFD(steps, length, h);
    u = est.getU();
end


function u = uXY(x, y)
    u = sin(x) + sin(y);
end