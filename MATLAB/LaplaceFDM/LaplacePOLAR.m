classdef LaplacePOLAR

    properties
        radius
        rStep
        tStep
        steps
        BCtype
    end

    methods

        % Constructor
        % Possibly not want complete circle or not at origin

        function lpcP = LaplacePOLAR(radius, steps, BCtype) 
            lpcP.radius = radius;
            lpcP.steps = steps;
            lpcP.BCtype = BCtype;
            lpcP.rStep = radius/steps(1);
            lpcP.tStep = (2*pi)/steps(2);
        end

        % ----------------------------------------

        % Method to get U values

        function U = runSim(lpcP)
            if lpcP.BCtype == 'd'
                A = lpcP.makeA();
                b = lpcP.makeb();
                x = linsolve(A,b);
                U = x;
            else
                error("Not made yet");
            end
        end

        % ----------------------------------------

        % Method to make A matrix to solve

        function A = makeA(lpcP)
            if lpcP.BCtype == 'd'
                A = zeros((lpcP.steps(1)+1)*lpcP.steps(2), (lpcP.steps(1)+1)*lpcP.steps(2));
                N = lpcP.steps(1)+1;
                M = lpcP.steps(2);
                NM = (lpcP.steps(1)+1)*lpcP.steps(2);
                for i = 1:NM
                    if mod(i, N) == 1 % Check This
                        A(i, i) = -2;
                        A(i, i+1) = 1;
                        A(i, i+N) = 1;
                    elseif i > 1 && i <= N
                        rs = (lpcP.rStep*i)^2;
                        A(i,i) = (-4*rs*(lpcP.tStep^2)-(4*(lpcP.rStep^2)))/(2*rs*(lpcP.rStep^2)*(lpcP.tStep^2)); % i == j
                        A(i,i+1) = (2*(lpcP.rStep*i)-lpcP.rStep)/(2*(lpcP.rStep*i)*rs); % i == j+1
                        A(i,i-1) = (2*(lpcP.rStep*i)+lpcP.rStep)/(2*(lpcP.rStep*i)*rs); %i == j-1
                        A(i,i+N) = 1/(((lpcP.rStep*i)^2)*((lpcP.tStep)^2)); % i == j+lpcP.steps(1)
                    elseif i > (M-1)*(N) && i <= NM
                        rs = (lpcP.rStep*i)^2;
                        A(i,i) = (-4*rs*(lpcP.tStep^2)-(4*(lpcP.rStep^2)))/(2*rs*(lpcP.rStep^2)*(lpcP.tStep^2)); % i == j
                        A(i,i+1) = (2*(lpcP.rStep*i)-lpcP.rStep)/(2*(lpcP.rStep*i)*rs); % i == j+1
                        A(i,i-1) = (2*(lpcP.rStep*i)+lpcP.rStep)/(2*(lpcP.rStep*i)*rs); %i == j-1
                        A(i,i-N) = 1/(((lpcP.rStep*i)^2)*((lpcP.tStep)^2)); % i == j-lpcP.steps(1)
                    else
                        rs = (lpcP.rStep*i)^2;
                        A(i,i) = (-4*rs*(lpcP.tStep^2)-(4*(lpcP.rStep^2)))/(2*rs*(lpcP.rStep^2)*(lpcP.tStep^2)); % i == j
                        A(i,i+1) = (2*(lpcP.rStep*i)-lpcP.rStep)/(2*(lpcP.rStep*i)*rs); % i == j+1
                        A(i,i-1) = (2*(lpcP.rStep*i)+lpcP.rStep)/(2*(lpcP.rStep*i)*rs); %i == j-1
                        A(i,i+N) = 1/(((lpcP.rStep*i)^2)*((lpcP.tStep)^2)); % i == j+lpcP.steps(1) 
                        A(i,i-N) = 1/(((lpcP.rStep*i)^2)*((lpcP.tStep)^2)); % i == j-lpcP.steps(1)
                    end
                end
            else
                error("Not made yet");
            end
        end

        % ----------------------------------------

        % Method to make b vector to solve

        function b = makeb(lpcP)
            if lpcP.BCtype == 'd'
                NM = (lpcP.steps(1)+1)*lpcP.steps(2);
                N = lpcP.steps(1)+1;
                M = lpcP.steps(2);
                b = zeros(NM, 1);
                cp = lpcP.centerPoint();
                for i = 1:NM % Check Coeffiecents
                    if i > 1 && i <= N
                        rs = (lpcP.rStep*i)^2;
                        p = (2*(lpcP.rStep*i)-lpcP.rStep)/(2*(lpcP.rStep*i)*rs);
                        b(i, 1) = -p*cp; 
                    elseif i > (M-1)*N && i <= NM
                        rs = (lpcP.rStep*i)^2;
                        p = (2*(lpcP.rStep*i)+lpcP.rStep)/(2*(lpcP.rStep*i)*rs);
                        b(i, 1) = -p*lpcP.gTheta(i*lpcP.tStep); 
                    end
                end
            else
                error("Not made yet");
            end
        end

        % ----------------------------------------

        % Method to take average of BC

        function ave = centerPoint(lpcP)
            ave = 0;
            for i = 1:lpcP.steps(2)
                ave = ave + lpcP.gTheta(i*lpcP.tStep);
            end
            ave = ave/lpcP.steps(2);
        end

        % ----------------------------------------

        % BC

        function val = gTheta(~, theta)
            val = sin(theta)+1;
        end

        % ----------------------------------------

    end
end