
bSteps = 79; % Amount of steps on the boundary -> only do odd... Explain in notes
tSteps = 500; % Amount of time steps
dr = (2*pi)/(bSteps+1); % Delta r
dtheta = (2*pi)/(bSteps+1); % Delta theta
k = 1e1; % Diffusivity Constant
dt = (k/2)*((dr^2)+(dtheta^2)) * .001; % So called CFL for Delta time

u = zeros(bSteps+2, bSteps+1, tSteps); % Initialize tensor u(r,theta,t)

for i = 1:bSteps+2 % Loop to traverse r
    for j = 1:bSteps+1 % Loop to traverse theta
        u(i,j,1) = IC((i-1)*dr, (j-1)*dtheta); % Enter Initial conditions to tensor
    end
end

for t = 2:tSteps % Loop to traverse time
    for i = bSteps+2:-1:1 % Loop to traverse r
        for j = 1:bSteps+1 % Loop to traverse theta
            ri = dr*(i-1); % Calculate ri to simplify calculations, subtract 1 to account for Matlab indexing

            % Constants: Refer to notes
            psi = ((-2*k*dt)*(((ri^2)*(dtheta^2))+(dr^2)))/((ri^2)*(dr^2)*(dtheta^2));  
            alpha = ((2*ri+dr)*k*dt)/(2*ri*(dr^2));
            beta = ((2*ri-dr)*k*dt)/(2*ri*(dr^2));
            phi = (k*dt)/((ri^2)*(dtheta^2));

            if i == 1 % if at the center -> u(0,theta)
                u(i,j,t) = (((-4*dt*k)/(dr^2))+1)*u(1,1,t-1) + ((dt*k)/(dr^2))*(u(2,1,t-1)+u(2,ceil(bSteps/4)+1,t-1)+u(2,ceil(bSteps/2)+1,t-1)+u(2,bSteps+1,t-1));
            elseif i == bSteps+2 && j == 1 % if at South-East edge -> Use given BC and north = south; 2nd and 5th term
                u(i, j, t) = ((psi+alpha+1)*u(i, j, t-1)) + (alpha*dr*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, bSteps+1, t-1));
            elseif i == bSteps+2 && j == bSteps+1 % if at North-East edge -> Use given BC and north = south; 2nd and 4th term
                u(i, j, t) = ((psi+alpha+1)*u(i, j, t-1)) + (alpha*dr*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, 1, t-1)) + (phi*u(i, j-1, t-1));
            elseif i == bSteps+2 % if your on East edge of r -> Use given BC function; look at 2nd term
                u(i, j, t) = ((psi+alpha+1)*u(i, j, t-1)) + (alpha*dr*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, j-1, t-1));
            elseif j == 1 % if your on the South edge of theta -> north edge = south edge; look at 5th term
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, bSteps+1, t-1));
            elseif j == bSteps+1 % if your on the North edge of theta -> north edge = south edge; look at 4th term
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, 1, t-1)) + (phi*u(i, j-1, t-1)); 
            else % if your in the middle of the grid
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, j-1, t-1));
            end
        end
    end
end

% Show tensor
% disp(u)

% Generate theta and r values based on the size of the matrix
[i, j, num_slices] = size(u);
theta = linspace(0, 2*pi, j);
r = linspace(0, 1, i);

% Convert polar coordinates to Cartesian coordinates for all slices
[Theta, R] = meshgrid(theta, r);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% Create a figure for plotting
figure;

%{

% Plot the slice at the beginning
subplot(1, 3, 1);
surf(X, Y, u(:,:,1));
xlabel('X');
ylabel('Y');
zlabel('U');
title('Slice at the Beginning');

% Plot the slice in the middle
middle_index = round(num_slices / 2);
subplot(1, 3, 2);
surf(X, Y, u(:,:,middle_index));
%zlim([-1,2]);
xlabel('X');
ylabel('Y');
zlabel('U');
title('Slice in the Middle');

% Plot the slice at the end
subplot(1, 3, 3);
surf(X, Y, u(:,:,end));
%zlim([-1,2]);
xlabel('X');
ylabel('Y');
zlabel('U');
title('Slice at the End');

%}

plt(X,Y,u,tSteps)


% Initial Condition function
function f = IC(r,theta)
f = cos(r); 
end

% Boundary Condition function
function g = BC(t,theta)
g = (cos(theta)/((t^2) + 10)); 
end

function plt(X, Y, TP, tIter)

for n=1:tIter
    surf(X, Y, TP(:, :, n));
    Mov(n)=getframe;
end

movie(Mov,1,1);
end




