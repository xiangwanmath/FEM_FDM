steps = 30;
dt = (2*pi)/steps;
dr = 1;

A = makeA(dr,dt, steps);
b = makeb(dr, steps);
U = jacobi(A,b,1e-5,100);

disp(U(1:steps+1));
%disp(A);
%disp(b);
%disp(U);



function A = makeA(dr,dt,steps)
n = steps+1;
A = zeros(n^2,n^2);
iC = 1;
jC = 1;
for i = 1:n
    for j = 1:n
        if i == 1 && j == 1
            A(iC:iC+n-1, jC:jC+n-1) = makeT0(steps);
        elseif i == j
            A(iC:iC+n-1, jC:jC+n-1) = makeT(i-1, dr, dt, steps);
        elseif i == j-1 && i ~= 1
            A(iC:iC+n-1, jC:jC+n-1) = makeIR(i-1, dr, dt, steps);
        elseif i == j+1 && i ~= 1
            A(iC:iC+n-1, jC:jC+n-1) = makeIL(i-1, dr, dt, steps);
        end
        jC = jC + n;
    end
    iC = iC + n;
    jC = 1;
end
end

function t0 = makeT0(steps)
t0 = eye(steps+1);
nV = -1*ones(steps,1);
t0 = t0 + diag(nV, 1);
t0(end,1) = -1;
end

function t = makeT(i, dr, dt, steps)
t = zeros(steps+1, steps+1);
iN = i * steps+1;
count = 1;
for j = iN+1:iN+steps+1
    ri = j * dr;
    psi = -2*((ri^2)*(dt^2) + dr^2)/((ri^2)*(dr^2)*(dt^2));
    t(count,count) = psi;
    phi = 1/((ri^2)*(dt^2));
    if count ~= 1
        t(count, count-1) = phi;
    end
    if count ~= steps+1
        t(count, count+1) = phi;
    end
    count = count + 1;
end
end

function ir = makeIR(i, dr, dt, steps)
ir = zeros(steps+1, steps+1);
iN = i * steps+1;
count = 1;
for j = iN+1:iN+steps+1
    ri = j * dr;
    alpha = (2*(ri^2) + ri*dr)/(2*(ri^2)*dt^2);
    ir(count,count) = alpha;
    count = count + 1;
end
end

function il = makeIL(i, dr, dt, steps)
il = zeros(steps+1, steps+1);
iN = i * steps+1;
count = 1;
for j = iN+1:iN+steps+1
    ri = j * dr;
    beta = (2*(ri^2) - ri*dr)/(2*(ri^2)*dt^2);
    il(count,count) = beta;
    count = count + 1;
end
end

function b = makeb(dr, steps)
b = zeros((steps+1)^2,1);
ri = dr * (steps-1); 
alpha = (2*(ri^2) + ri*dr)/(2*(ri^2)*(dr^2));
for i = 0:steps
    b(((steps+1)*steps)+i,1) = -alpha*f(i*dr);
end
end

function u = jacobi(A,b,tol,maxIts)
end

function f = f(theta)
f = 4*sin(theta);
end




