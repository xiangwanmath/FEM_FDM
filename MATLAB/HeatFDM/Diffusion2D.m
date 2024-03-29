classdef Diffusion2D

    properties
        dt
        dx
        dy
        xSteps
        ySteps
        tIter
        BCtype
        BCoord
        k
    end

    methods

        % Constructor

        function dif = Diffusion2D(tStep, tIter, xSteps, ySteps, k, boundCoord, boundCondsType)
            dif.dt = tStep;
            dif.xSteps = xSteps;
            dif.ySteps = ySteps;
            dif.tIter = tIter;
            dif.k = k;
            dif.BCoord = boundCoord;
            dif.BCtype = boundCondsType;
            dif.dx = abs(boundCoord(1,2)-boundCoord(1,1))/xSteps; 
            dif.dy = abs(boundCoord(2,2)-boundCoord(2,1))/ySteps;
        end

        % -------------------------------

        % Function to run the model
        % Optimization:
        % Use Gauss-Siedel method

        function u = runSim(dif)
            if strcmp(dif.BCtype(1),'d') && strcmp(dif.BCtype(2), 'd') && strcmp(dif.BCtype(3), 'd') && strcmp(dif.BCtype(4), 'd')
                u = zeros(dif.xSteps, dif.ySteps, dif.tIter+1);
                A = dif.makeA();
                t = 0;
                u(:, :, 1) = dif.makeIC();
                for i = 1:dif.tIter
                    t = t + dif.dt;
                    b = dif.makeb(t, u(:,:, i));
                    uTempV = linsolve(A,b); 
                    uTempM = zeros(dif.xSteps,dif.ySteps);
                    count = 1;
                    for j = 1:dif.xSteps
                        for l = 1:dif.ySteps
                            uTempM(j, l) = uTempV(count);
                            count = count + 1;
                        end
                    end
                    u(:,:,i+1) = uTempM; 
                end
            else
                error("Not made yet");
            end
        end

        % -------------------------------

        % Function to create A in Ax=b:
        % Errors:
        % 0s on diag?
        % Optimization:
        % Possibly use sparse function

        function A = makeA(dif)
            xCon = (dif.k*dif.dt)/(dif.dx^2);
            yCon = (dif.k*dif.dt)/(dif.dy^2);
            if dif.BCtype(1) == 'd' && dif.BCtype(2) == 'd' && dif.BCtype(3) == 'd' && dif.BCtype(4) == 'd'
                mid = repmat((1+(2*xCon)+(2*yCon)), dif.xSteps^2, 1);
                sX = repmat(-xCon, (dif.xSteps^2)-1,1);
                sY = repmat(-yCon, (dif.xSteps^2)-dif.xSteps,1);
                A = diag(mid) + diag(sX,1) + diag(sX,-1) + diag(sY,dif.xSteps) + diag(sY,-dif.xSteps);
            else
                error("Not made yet");
            end
        end

        % -------------------------------

        % Function to create b in Ax=b at some time t:

        function b = makeb(dif, t, prev)
            if dif.BCtype(1) == 'd' && dif.BCtype(2) == 'd' && dif.BCtype(3) == 'd' && dif.BCtype(4) == 'd'
                b = zeros(dif.xSteps*dif.ySteps, 1);
                count = 1;
                xCon = (dif.k*dif.dt)/(dif.dx^2);
                yCon = (dif.k*dif.dt)/(dif.dy^2);
                for i = 1:dif.xSteps
                    for j = 1:dif.ySteps

                        if i == 1 && j == 1
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gW(dif.dy*j, t) + xCon*dif.gS(dif.dx*i, t);
                        elseif i == 1 && j == dif.ySteps
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gW(dif.dy*j, t) + xCon*dif.gN(dif.dx*i, t);
                        elseif i == dif.xSteps && j == 1
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gE(dif.dy*j, t) + xCon*dif.gS(dif.dx*i, t);
                        elseif i == dif.xSteps && j == dif.ySteps
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gE(dif.dy*j, t) + xCon*dif.gN(dif.dx*i, t);
                        elseif i == 1
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gW(dif.dy*j, t);
                        elseif j == 1
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + xCon*dif.gS(dif.dx*i, t);
                        elseif i == dif.xSteps
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + yCon*dif.gE(dif.dy*j, t);
                        elseif j == dif.ySteps
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt + xCon*dif.gN(dif.dx*i, t);
                        else
                            b(count, 1) = prev(i,j) + dif.fXYT((i*dif.dx)+dif.BCoord(1,1),(j*dif.dy)+dif.BCoord(2,1),t)*dif.dt;
                        end

                        count = count + 1;
                    end
                end
            else
                error("Not made yet");
            end
        end

        % -------------------------------

        % Function for Gauss-Siedel 
        % Check Laplace

        function m = gaussSiedel(~, A, b)
        end

        % -------------------------------

        % Function to tell if a sequence of vectors has converged

        function bool = hasConverged(~, vNext, vCurr, percent)
            numer = norm(vNext-vCurr);
            denom = norm(vCurr);
            error = (numer/denom)*100;
            if error <= percent
                bool = true;
            else 
                bool = false;
            end
        end

        % -------------------------------

        % Function to get matrix of IC values

        function IC = makeIC(dif)
            IC = zeros(dif.xSteps, dif.ySteps); 
            for i = 1:dif.xSteps
                for j = 1:dif.ySteps
                    IC(i, j) = dif.uXY0(((i*dif.dx)+dif.BCoord(1,1)),((j*dif.dy)+dif.BCoord(2,1)));
                end
            end
        end

        % -------------------------------

        % Functions to create mesh grids 

        function x = getXGrid(dif)
            x = zeros(dif.xSteps, dif.ySteps);
            for i = 1:dif.xSteps
                x(:, i) = (i*dif.dx) + dif.BCoord(1,1);
            end
        end

        function y = getYGrid(dif)
            y = zeros(dif.xSteps, dif.ySteps);
            for i = 1:dif.ySteps
                y(i, :) = (i*dif.dy) + dif.BCoord(2,1);
            end
        end

        % -------------------------------

        % Initial Condition Function: u(x,y,0)

        function uXY0 = uXY0(dif, x, y) 
            uXY0 = (2*dif.k + 1)*cos(x)*cos(y);
        end

        % -------------------------------

        % Source Function: ut = uxx + f(x,y,t)
        % Assuming the solution is u(x,y,t) = e^-t * cos(x) * cos(y)

        function fXYT = fXYT(dif, x, y, t)
            fXYT = (2*dif.k - 1)*exp(-t)*cos(x)*cos(y);
        end

        % -------------------------------

        % Boundary Condition Functions:
        % N (North) -> BC @ YL
        % E (East) -> BC @ XL
        % S (South) -> BC @ Y0 
        % W (West) -> BC @ X0 

        function gN = gN(~, x, t)
            gN = 0;
        end

        function gE = gE(~, y, t)
            gE = 0;
        end

        function gS = gS(~, x, t)
            gS = 0;
        end

        function gW = gW(~, y, t)
            gW = 0;
        end

        % --------------------------------
    end
end