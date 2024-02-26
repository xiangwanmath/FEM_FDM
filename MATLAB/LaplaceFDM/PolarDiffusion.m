bSteps = 10;
tSteps = 10;
dr = (2*pi)/(bSteps+1);
dtheta = (2*pi)/(bSteps+1);
k = 1e-5;
dt = (k/2)*((dr^2)+(dtheta^2));

u = zeros(bSteps+1, bSteps+1, tSteps);

for i = 1:bSteps+1
    for j = 1:bSteps+1
        u(i,j,1) = IC((i-1)*dr, (j-1)*dtheta);
    end
end

for t = 2:tSteps
    for i = bSteps+1:-1:1
        for j = 1:bSteps+1
            ri = dr*(i-1);
            if i == 1
                ri = dr*i;
            end
            psi = ((-2*k*dt)*(((ri^2)*(dtheta^2))+(dr^2)))/((ri^2)*(dr^2)*(dtheta^2));  
            alpha = ((2*ri+dr)*k*dt)/(2*ri*(dr^2));
            beta = ((2*ri-dr)*k*dt)/(2*ri*(dr^2));
            phi = (k*dt)/((ri^2)*(dtheta^2));
            if i == 1
                jUp = j+1;
                jDown = j-1;
                if j == 1
                    jDown = bSteps+1;
                elseif j == bSteps+1
                    jUp = 1;
                end
                u(i, j, t) = (1/beta)*(((1/dt)*(u(i+1,j,t)-u(i+1,j,t-1)))-psi*u(i+1,j,t)-alpha*u(i+2,j,t)-phi*u(i+1,jUp,t)-phi*u(i+1,jDown,t));
                %disp(u(i,j,t)); % Diverging...
            elseif i == bSteps+1 && j == 1
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, bSteps+1, t-1));
                %disp(u(i,j,t)); % What
            elseif i == bSteps+1 && j == bSteps+1
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, 1, t-1)) + (phi*u(i, j-1, t-1));
                %disp(u(i,j,t));
            elseif i == bSteps+1
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*BC(dt*t,dtheta*j)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, j-1, t-1));
                %disp(u(i,j,t));
            elseif j == 1
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, bSteps+1, t-1));
                %disp(u(i,j,t));
            elseif j == bSteps+1
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, 1, t-1)) + (phi*u(i, bSteps+1, t-1));
            else
                u(i, j, t) = ((psi+1)*u(i, j, t-1)) + (alpha*u(i+1, j, t-1)) + (beta*u(i-1, j, t-1)) + (phi*u(i, j+1, t-1)) + (phi*u(i, j-1, t-1));
            end
        end
    end
end

disp(u)

function f = IC(r,theta)
f = sin(r)*sin(r); 
end

function g = BC(t,theta)
g = 0; 
end



