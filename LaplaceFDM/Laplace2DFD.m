classdef Laplace2DFD

    properties
        steps = []
        lengths = []
        h
    end

    methods

        function lpc = Laplace2DFD(steps, lengths, h)
            lpc.steps = steps;
            lpc.lengths = lengths;
            lpc.h = h;
        end

        function U = getU(lpc)
            A = lpc.makeA(); 
            f = lpc.makef();
            U = linsolve(A, f);
        end

        function uMat = runSimNeu(lpc)
            A = lpc.makeANeu(); 
            f = lpc.makefNeu();
            U = linsolve(A, f);
            uMat = zeros(lpc.steps(1)+2, lpc.steps(1)+2);
            for i = 0:lpc.steps(1)+1
                uMat(i+1, 1:lpc.steps(1)+2) = U((i*(lpc.steps(1)+2))+1:(i*(lpc.steps(1)+2))+lpc.steps(1)+2);
            end
        end

        function uMat = runSim(lpc)
            A = lpc.makeA(); 
            f = lpc.makef();
            U = linsolve(A, f);
            uMat = zeros(lpc.steps(1), lpc.steps(1));
            for i = 0:lpc.steps(1)+1
                uMat(i+1, 1:lpc.steps(1)+1) = U(i*lpc.steps(1)+1:i*lpc.steps(1)+lpc.steps(1)); 
            end
        end


        function A = makeA(lpc) 
            mid = repmat(-4/lpc.h^2, lpc.steps(1)^2, 1);
            id = repmat(1/lpc.h^2, (lpc.steps(1)^2)-lpc.steps(1), 1);
            onesV = repmat(1/lpc.h^2, lpc.steps(1)-1, 1);
            onesV = [onesV; 0];
            ones = repmat(onesV, lpc.steps(1), 1);
            ones(lpc.steps(1)^2) = [];
            A = diag(mid) + diag(ones, -1) + diag(ones, 1) + diag(id, lpc.steps(1)) + diag(id, -lpc.steps(1));
        end

        function A = makeANeu(lpc)
            mid = repmat(-4/lpc.h^2, (lpc.steps(1)+2)^2, 1);
            id = repmat(1/lpc.h^2, ((lpc.steps(1)+2)^2)-(lpc.steps(1)+2), 1);
            onesV = repmat(1/lpc.h^2, lpc.steps(1)+1, 1);
            onesV = [onesV; 0];
            onesV(1,1) = 2/lpc.h^2;
            onesV(lpc.steps(1)+1,1) = 2/lpc.h^2;
            ones = repmat(onesV, lpc.steps(1)+2, 1);
            ones((lpc.steps(1)+2)^2) = [];
            A = diag(mid) + diag(ones, -1) + diag(ones, 1) + diag(id, lpc.steps(1)+2) + diag(id, -(lpc.steps(1)+2));
        end

        function f = makef(lpc)
            f = zeros(lpc.steps(1)^2, 1);
            hP = 1/(lpc.h^2);
            count = 1;
            for i = 1:lpc.steps(1)
                for j = 1:lpc.steps(1)
                    if i == 1
                        if j == 1
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCX0(i*lpc.h)+lpc.uBCY0(j*lpc.h)));
                        elseif j == lpc.steps(1)
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCXL(i*lpc.h)+lpc.uBCY0(j*lpc.h)));
                        else
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCY0(j*lpc.h)));
                        end
                    elseif i == lpc.steps(1)
                        if j == 1
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCX0(i*lpc.h)+lpc.uBCYL(j*lpc.h)));
                        elseif j == lpc.steps(1)
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCXL(i*lpc.h)+lpc.uBCYL(j*lpc.h)));
                        else
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCYL(j*lpc.h)));
                        end
                    else
                        if j == 1
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCX0(i*lpc.h)));
                        elseif j == lpc.steps(1)
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - (hP*(lpc.uBCXL(i*lpc.h)));
                        else
                            temp = lpc.fXY(i*lpc.h, j*lpc.h);
                        end
                    end 
                    f(count, 1) = temp;
                    count = count + 1;
                end
            end
        end

        function f = makefNeu(lpc)
            f = zeros((lpc.steps(1)+2)^2, 1);
            hP = 2/lpc.h;
            count = 1;
            for i = 0:lpc.steps(1)+1
                for j = 0:lpc.steps(1)+1
                    if i == 0
                        if j == 0
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) + hP*lpc.gW(i*lpc.h) + hP*lpc.gS(j*lpc.h);
                        elseif j == lpc.steps(1)+1
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) + hP*lpc.gS(i*lpc.h) - hP*lpc.gN(j*lpc.h);
                        else
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) + hP*lpc.gW(j*lpc.h);
                        end
                    elseif i == lpc.steps(1)+1
                        if j == 0
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - hP*lpc.gE(i*lpc.h) - hP*lpc.gS(j*lpc.h);
                        elseif j == lpc.steps+1
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - hP*lpc.gE(i*lpc.h) - hP*lpc.gN(j*lpc.h);
                        else
                            temp = lpc.fXY(i*lpc.h, j*lpc.h) - hP*lpc.gE(j*lpc.h);
                        end
                    elseif j == 0
                        temp = lpc.fXY(i*lpc.h, j*lpc.h) + hP*lpc.gS(i*lpc.h);
                    elseif j == lpc.steps(1)+1
                         temp = lpc.fXY(i*lpc.h, j*lpc.h) - hP*lpc.gN(i*lpc.h);
                    else
                        temp = lpc.fXY(i*lpc.h, j*lpc.h);
                    end
                    f(count, 1) = temp;
                    count = count + 1;
                end
            end
        end

        function uBCX0 = uBCX0(~, val)
            uBCX0 = sin(val);
        end

        function uBCXL = uBCXL(~, val)
            uBCXL = sin(val);
        end

        function uBCY0 = uBCY0(~, val)
            uBCY0 = sin(val);
        end

        function uBCYL = uBCYL(~, val)
            uBCYL = sin(val);
        end

        function uxBCX0 = gW(~, val)
            uxBCX0 = cos(val);
        end

        function uxBCXL = gE(~, val)
            uxBCXL = sin(val);
        end

        function uyBCY0 = gS(~, val)
            uyBCY0 = sin(val);
        end

        function uyBCYL = gN(~, val)
            uyBCYL = sin(val);
        end

        function fXY = fXY(~, x, y)
            fXY = -sin(x)-sin(y) + tan(y^2);
        end

    end
end