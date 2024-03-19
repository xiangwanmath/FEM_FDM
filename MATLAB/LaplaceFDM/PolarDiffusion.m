
bSteps = 15; % Amount of steps on the boundary -> only do odd... and -1 must must divisible by 4Explain in notes
tSteps = 100; % Amount of time steps
dr = (pi)/(bSteps+1); % Delta r
dtheta = (2*pi)/(bSteps+1); % Delta theta
k = 1e1; % Diffusivity Constant
dt = (k/2)*((dr^2)+(dtheta^2)) * .001; % So called CFL for Delta time

u = zeros(bSteps+1, bSteps+1, tSteps); % Initialize tensor u(x,y,t)

for i = 1:bSteps+1 % Loop to traverse r
    for j = 1:bSteps+1 % Loop to traverse theta
        u(i,j,1) = IC((i-1)*dr, (j-1)*dtheta); % Enter Initial conditions to tensor
    end
end

for t = 2:tSteps % Loop to traverse time
    for i = bSteps+1:-1:1 % Loop to traverse r
        for j = 1:bSteps+1 % Loop to traverse theta
            ri = dr*(i-1); % Calculate ri to simplify calculations, subtract 1 to account for Matlab indexing

            % Constants: Refer to notes
            psi = ((-2*k*dt)*(((ri^2)*(dtheta^2))+(dr^2)))/((ri^2)*(dr^2)*(dtheta^2));  
            alpha = ((2*ri+dr)*k*dt)/(2*ri*(dr^2));
            beta = ((2*ri-dr)*k*dt)/(2*ri*(dr^2));
            phi = (k*dt)/((ri^2)*(dtheta^2));

            if i == 1 % if at the center -> u(0,theta)
                u(i,j,t) = ((-4*dt*k)+1)*u(1,1,t-1) + (-4*dt*k)*(u(2,1,t-1)+u(2,ceil(bSteps/4)+1,t-1)+u(2,ceil(bSteps/2)+1,t-1)+u(2,bSteps+1,t-1));
            elseif i == bSteps+1 && j == 1 % if at South-East edge -> Use given BC and north = south; 2nd and 5th term
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, bSteps+1, t-1));
            elseif i == bSteps+1 && j == bSteps+1 % if at North-East edge -> Use given BC and north = south; 2nd and 4th term
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, 1, t-1)) + (phi*u(i, j-1, t-1));
            elseif i == bSteps+1 % if your on East edge of r -> Use given BC function; look at 2nd term
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, j-1, t-1));
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
disp(u)

% Assuming you have your matrix of u values named 'u'

% Generate theta and r values based on the size of the matrix
[num_rows, num_cols, num_slices] = size(u);
theta = linspace(0, 2*pi, num_rows);
r = linspace(0, 1, num_cols);

% Convert polar coordinates to Cartesian coordinates for all slices
[Theta, R] = meshgrid(theta, r);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% Create a figure for plotting
figure;

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
xlabel('X');
ylabel('Y');
zlabel('U');
title('Slice in the Middle');

% Plot the slice at the end
subplot(1, 3, 3);
surf(X, Y, u(:,:,end));
xlabel('X');
ylabel('Y');
zlabel('U');
title('Slice at the End');


% Initial Condition function
function f = IC(r,theta)
f = sin(r); 
end

% Boundary Condition function
function g = BC(t,theta)
g = 0; 
end




