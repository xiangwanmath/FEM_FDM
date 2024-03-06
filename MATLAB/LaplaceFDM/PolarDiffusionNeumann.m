
bSteps = 11; % Amount of steps on the boundary -> only do odd... Explain in notes
tSteps = 10; % Amount of time steps
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
            % ^ think about this more

            if i == 1 % 1 is 0... but we are using r1 to solve for center
                ri = dr*i;
            end

            % Constants: Refer to notes
            psi = ((-2*k*dt)*(((ri^2)*(dtheta^2))+(dr^2)))/((ri^2)*(dr^2)*(dtheta^2));  
            alpha = ((2*ri+dr)*k*dt)/(2*ri*(dr^2));
            beta = ((2*ri-dr)*k*dt)/(2*ri*(dr^2));
            phi = (k*dt)/((ri^2)*(dtheta^2));

            if i == 1 % if at the center -> u(0,theta)
                tN = j+((bSteps+1)/2);

                if tN > bSteps+1
                    tN = mod(tN, bSteps+1);
                end

                % Do this conditional so you don't cause indexing errors at
                % corners
                jDown = j-1;
                jUp = j+1;
                if j == 1
                    jDown = bSteps+1;
                elseif j == bSteps+1
                    jUp = 1;
                end
                %--------------------------------------------------------

                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i+1, tN, t-1)) + (phi*u(i, jUp, t-1)) + (phi*u(i, jDown, t-1));
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
disp(u)


% time_point = 3;
% disp((u(2,1,time_point) - u(2,1,time_point - 1))/dt)
% u(2,1,time_point - 1)
% u(2,1,time_point)


% Initial Condition function
function f = IC(r,theta)
f = sin(r)*sin(r); 
end

% Boundary Condition function
function g = BC(t,theta)
g = 0; 
end


