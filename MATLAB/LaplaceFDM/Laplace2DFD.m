classdef Laplace2DFD

    properties
        steps = []
        BCoords = []
        BCtype = []
        h
        alpha
    end

    methods

        % Constructor

        function lpc = Laplace2DFD(steps, BCoords, BCtype)
            lpc.steps = steps;
            lpc.BCoords = BCoords;
            lpc.h = abs(BCoords(1,2) - BCoords(1,1))/steps(1);
            lpc.BCtype = BCtype;
            lpc.alpha = 6;
        end

        % ----------------------------

        % Run Sim

        function uMat = runSim(lpc)
            A = lpc.makeA();
            b = lpc.makeb();
            %U = linsolve(A, b);
            len = size(b);
            x0 = zeros(len(1),1);
            x0(1, 1) = 1;
            U = lpc.jacobi(A, b, x0, 1e-6, 100);
            if lpc.BCtype(1) == "d" && lpc.BCtype(2) == "d" && lpc.BCtype(3) == "d" && lpc.BCtype(4) == "d"
                uMat = zeros(lpc.steps(1), lpc.steps(1));
                count = 1;
                for i = 1:lpc.steps(1)
                    for j = 1:lpc.steps(2)
                        uMat(i, j) = U(count);
                        count = count + 1;
                    end
                end
            elseif (lpc.BCtype(1) == "n" && lpc.BCtype(2) == "n" && lpc.BCtype(3) == "n" && lpc.BCtype(4) == "n") || (lpc.BCtype(1) == "r" && lpc.BCtype(2) == "r" && lpc.BCtype(3) == "r" && lpc.BCtype(4) == "r")
                uMat = zeros(lpc.steps(1)+2, lpc.steps(1)+2);
                count = 1;
                for i = 1:lpc.steps(1)+2
                    for j = 1:lpc.steps(2)+2
                        uMat(i,  j) = U(count);
                        count = count + 1;
                    end
                end
            else
                error("Not Made Yet")
            end
        end

        % ----------------------------

        % Make A

        function A = makeA(lpc)
            if lpc.BCtype(1) == "d" && lpc.BCtype(2) == "d" && lpc.BCtype(3) == "d" && lpc.BCtype(4) == "d"
                mid = repmat(-4/lpc.h^2, lpc.steps(1)^2, 1);
                id = repmat(1/lpc.h^2, (lpc.steps(1)^2)-lpc.steps(1), 1);
                onesV = repmat(1/lpc.h^2, lpc.steps(1)-1, 1);
                onesV = [onesV; 0];
                ones = repmat(onesV, lpc.steps(1), 1);
                ones(lpc.steps(1)^2) = [];
                A = diag(mid) + diag(ones, -1) + diag(ones, 1) + diag(id, lpc.steps(1)) + diag(id, -lpc.steps(1));
            elseif lpc.BCtype(1) == "n" && lpc.BCtype(2) == "n" && lpc.BCtype(3) == "n" && lpc.BCtype(4) == "n"
                mid = repmat(-4/lpc.h^2, (lpc.steps(1)+2)^2, 1);
                id = repmat(1/lpc.h^2, ((lpc.steps(1)+2)^2)-(lpc.steps(1)+2), 1);
                onesV = repmat(1/lpc.h^2, lpc.steps(1)+1, 1);
                onesV = [onesV; 0];
                onesV(1,1) = 2/lpc.h^2;
                onesV(lpc.steps(1)+1,1) = 2/lpc.h^2;
                ones = repmat(onesV, lpc.steps(1)+2, 1);
                ones((lpc.steps(1)+2)^2) = [];
                A = diag(mid) + diag(ones, -1) + diag(ones, 1) + diag(id, lpc.steps(1)+2) + diag(id, -(lpc.steps(1)+2));
            elseif lpc.BCtype(1) == "r" && lpc.BCtype(2) == "r" && lpc.BCtype(3) == "r" && lpc.BCtype(4) == "r"
                mid = repmat((2-lpc.alpha)/lpc.h^2, (lpc.steps(1)+2)^2, 1);
                id = repmat(1/lpc.h^2, ((lpc.steps(1)+2)^2)-(lpc.steps(1)+2), 1);
                onesV = repmat(1/lpc.h^2, lpc.steps(1)+1, 1);
                onesV = [onesV; 0];
                onesV(1,1) = 2/lpc.h^2;
                onesV(lpc.steps(1)+1,1) = 2/lpc.h^2;
                ones = repmat(onesV, lpc.steps(1)+2, 1);
                ones((lpc.steps(1)+2)^2) = [];
                A = diag(mid) + diag(ones, -1) + diag(ones, 1) + diag(id, lpc.steps(1)+2) + diag(id, -(lpc.steps(1)+2));
            else
                error("Not Made Yet")
            end
        end

        % ----------------------------

        % Make b
        % Fix Points
        % Approximating cos(3.1416) = -1 ??

        function b = makeb(lpc)
            if lpc.BCtype(1) == "d" && lpc.BCtype(2) == "d" && lpc.BCtype(3) == "d" && lpc.BCtype(4) == "d"
                b = zeros(lpc.steps(1)^2, 1);
                hP = 1/(lpc.h^2);
                count = 1;
                for i = 1:lpc.steps(1)
                    for j = 1:lpc.steps(1)

                        if i == 1 && j == 1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gW((i*lpc.h)+lpc.BCoords(1,1))+lpc.gS((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif i == 1 && j == lpc.steps(1)
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gE((i*lpc.h)+lpc.BCoords(1,1))+lpc.gS((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif i == lpc.steps(1) && j == 1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gW((i*lpc.h)+lpc.BCoords(1,1))+lpc.gN((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif i == lpc.steps(1) && j == lpc.steps(1)
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gE((i*lpc.h)+lpc.BCoords(1,1))+lpc.gN((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif i == 1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gS((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif j == 1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gW((i*lpc.h)+lpc.BCoords(1,1))));
                        elseif i == lpc.steps(1)
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gN((j*lpc.h)+lpc.BCoords(2,1))));
                        elseif j == lpc.steps(1)
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - (hP*(lpc.gE((i*lpc.h)+lpc.BCoords(1,1))));
                        else
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1));
                        end

                        count = count + 1;
                    end
                end
            elseif lpc.BCtype(1) == "n" && lpc.BCtype(2) == "n" && lpc.BCtype(3) == "n" && lpc.BCtype(4) == "n"
                b = zeros((lpc.steps(1)+2)^2, 1);
                hP = 2/lpc.h;
                count = 1;
                for i = 0:lpc.steps(1)+1
                    for j = 0:lpc.steps(1)+1 

                        if i == 0 && j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gW((i*lpc.h)+lpc.BCoords(1,1)) + hP*lpc.gS((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == 0 && j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gS((i*lpc.h)+lpc.BCoords(1,1)) - hP*lpc.gN((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gW((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1 && j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((i*lpc.h)+lpc.BCoords(1,1)) + hP*lpc.gS((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1 && j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((i*lpc.h)+lpc.BCoords(1,1)) - hP*lpc.gN((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((j*lpc.h)+lpc.BCoords(2,1));
                        elseif j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gS((i*lpc.h)+lpc.BCoords(1,1));
                        elseif j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gN((i*lpc.h)+lpc.BCoords(1,1));
                        else
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1));
                        end
                        
                        count = count + 1;
                    end
                end
            elseif lpc.BCtype(1) == "r" && lpc.BCtype(2) == "r" && lpc.BCtype(3) == "r" && lpc.BCtype(4) == "r"
                b = zeros((lpc.steps(1)+2)^2, 1);
                hP = 2/lpc.h;
                count = 1;
                for i = 0:lpc.steps(1)+1
                    for j = 0:lpc.steps(1)+1 

                        if i == 0 && j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gW((i*lpc.h)+lpc.BCoords(1,1)) + hP*lpc.gS((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == 0 && j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gS((i*lpc.h)+lpc.BCoords(1,1)) - hP*lpc.gN((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gW((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1 && j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((i*lpc.h)+lpc.BCoords(1,1)) + hP*lpc.gS((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1 && j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((i*lpc.h)+lpc.BCoords(1,1)) - hP*lpc.gN((j*lpc.h)+lpc.BCoords(2,1));
                        elseif i == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gE((j*lpc.h)+lpc.BCoords(2,1));
                        elseif j == 0
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) + hP*lpc.gS((i*lpc.h)+lpc.BCoords(1,1));
                        elseif j == lpc.steps(1)+1
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1)) - hP*lpc.gN((i*lpc.h)+lpc.BCoords(1,1));
                        else
                            b(count, 1) = lpc.fXY((i*lpc.h)+lpc.BCoords(1,1), (j*lpc.h)+lpc.BCoords(2,1));
                        end
                        
                        count = count + 1;
                    end
                end
            else
                error("Not Made Yet");
            end
        end

        % ----------------------------

        % BC

        function gW = gW(~, val)
            gW = 0;
        end

        function gE = gE(~, val)
            gE = 0;
        end

        function gS = gS(~, val)
            gS = 0;
        end

        function gN = gN(~, val)
            gN = 0;
        end

        % ----------------------------

        % Source

        function fXY = fXY(~, x, y)
            fXY = -2*cos(x)*cos(y);
        end

        % ----------------------------

        function x = jacobi(~, A, b, x0, tol, max_iter)
            [~,n] = size(A);

            x = x0;
            iter = 0;
            converged = false;

            while ~converged && iter < max_iter
                x_old = x;

                for i = 1:n
                    sigma = A(i, :) * x_old - A(i, i) * x_old(i);
                    x(i) = (b(i) - sigma) / A(i, i);
                end

                if norm(x - x_old, inf) < tol
                    converged = true;
                end

                iter = iter + 1;
            end
        end

        % ----------------------------

    end
end