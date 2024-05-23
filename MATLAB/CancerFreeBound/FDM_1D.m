classdef FDM_1D
    
    properties

        tSteps
        rSteps
        
        C
        T
        W
        
        R
        t
        
        dr
        dt
        
        gammaC
        gammaT
        gammaW
        
        lambdaC
        lambdaT
        
        Wh
        W0
        K
        nu
        
        dT
        epsilon
        T0
        
        dWC
        dW
    end

%{
To do's:
Fill in BC/IC
Fill in Constants
Fill in Aux functions
Check all indexing
Test non-adaptable dt
Test adaptable dt
%}
    
    methods
        function cfb = FDM_1D(tSteps, rSteps, initR)

            cfb.tSteps = tSteps;
            cfb.rSteps = rSteps;
            
            cfb.tSteps = tSteps;
            cfb.rSteps = rSteps; % includes boundaries

            cfb.C = zeros(rSteps, tSteps);
            cfb.T = zeros(rSteps, tSteps);
            cfb.W = zeros(rSteps, tSteps);

            cfb.R = zeros(1, tSteps);
            cfb.R(1) = initR;
            cfb.t = zeros(1, tSteps);
            cfb.t(1) = 0;

            cfb.dr = zeros(1, tSteps);
            cfb.dr(1) = cfb.R(1)/(cfb.rSteps-1); % check this
            cfb.dt = zeros(1, rSteps);
            cfb.dt(1) = 0.1; % fix later

            cfb.gammaC = 8.64e-7;
            cfb.gammaT = 8.64e-7;
            cfb.gammaW = 2;

            cfb.lambdaC = 2.24;
            cfb.lambdaT = 4e-3;

            cfb.Wh = 2.5e-6;
            cfb.W0 = 4e-6;
            cfb.K = .4;
            cfb.nu = 100000; % can change

            cfb.dT = 0.18;
            cfb.epsilon = 0.1; % have to find
            cfb.T0 = 1e-3;

            cfb.dWC = 2.08;
            cfb.dW = 1.04;
        end

        function [t, R] = run_euler(cfb)

            for i = 1:cfb.rSteps
                cfb.C(i,1) = cfb.C_IC(cfb.dr(1)*(i-1));
                cfb.T(i,1) = cfb.T_IC(cfb.dr(1)*(i-1));
                cfb.W(i,1) = cfb.W_IC(cfb.dr(1)*(i-1));
            end

            cfb.C'
            cfb.T'
            cfb.W'

            for j = 2 : cfb.tSteps

                left = 0;
                right = 0;
                % quad = 0;
                disp("next time step")
                for i = 2:cfb.rSteps

                    sourceC = (cfb.lambdaC/(cfb.W0-cfb.Wh))*max(0,cfb.W(i,j-1)-cfb.Wh)*cfb.C(i,j-1)*(cfb.C(i,j-1)-(cfb.C(i,j-1)/cfb.K))-(cfb.nu*cfb.T(i,j-1)*cfb.C(i,j-1));
                    sourceT = cfb.lambdaT*cfb.C(i,j-1)-(cfb.dT*(cfb.T(i,j-1)-(cfb.epsilon*cfb.T0)));
                    sourceW = -(cfb.dWC*cfb.C(i,j-1)*cfb.W(i,j-1)) - (cfb.dW*cfb.W(i,j-1));

                    % disp(sourceC);
                    
                    alpha = (2*(i-1)+1)/(2*(i-1)*cfb.dr(j-1)*cfb.dr(j-1));
                    psi = (2*(i-1)-1)/(2*(i-1)*cfb.dr(j-1)*cfb.dr(j-1));
                    phi = (-2)/(cfb.dr(j-1)*cfb.dr(j-1));
                
                    if i == 1
                        cfb.C(i,j) = cfb.dt(j-1)*(cfb.gammaC*(((alpha+psi)*cfb.C(i+1,j-1)) + (phi*cfb.C(i,j-1))) + sourceC) + cfb.C(i, j-1);
                        cfb.T(i,j) = cfb.dt(j-1)*(cfb.gammaT*(((alpha+psi)*cfb.T(i+1,j-1)) + (phi*cfb.T(i,j-1))) + sourceT) + cfb.T(i, j-1);
                        cfb.W(i,j) = cfb.dt(j-1)*(cfb.gammaW*(((alpha+psi)*cfb.W(i+1,j-1)) + (phi*cfb.W(i,j-1))) + sourceW) + cfb.W(i, j-1);
                    elseif i == cfb.rSteps
                        alphaP = 1;
                        cfb.C(i,j) = cfb.C_BC(cfb.C(i,j-1), alphaP, cfb.dr(j-1));
                        cfb.T(i,j) = cfb.T_BC(cfb.T(i,j-1), cfb.T0, cfb.dr(j-1));
                        cfb.W(i,j) = cfb.W_BC(cfb.t(j-1));
                    else
                        cfb.C(i,j) = cfb.dt(j-1)*(cfb.gammaC*((alpha*cfb.C(i+1,j-1)) + (phi*cfb.C(i,j-1)) + (psi*cfb.C(i-1,j-1))) + sourceC) + cfb.C(i, j-1); 
                        cfb.T(i,j) = cfb.dt(j-1)*(cfb.gammaT*((alpha*cfb.T(i+1,j-1)) + (phi*cfb.T(i,j-1)) + (psi*cfb.T(i-1,j-1))) + sourceT) + cfb.T(i, j-1);
                        cfb.W(i,j) = cfb.dt(j-1)*(cfb.gammaW*((alpha*cfb.W(i+1,j-1)) + (phi*cfb.W(i,j-1)) + (psi*cfb.W(i-1,j-1))) + sourceW) + cfb.W(i, j-1);
                    end

                    
                    if i == 1
                        left = left + (sourceC*(cfb.dr(j-1)*(i-1))*(cfb.dr(j-1)*(i-1)));
                    elseif i == cfb.rSteps
                        right = right + (sourceC*(cfb.dr(j-1)*(i-1))*(cfb.dr(j-1)*(i-1)));
                    else
                        left = left + sourceC;
                        right = right + (sourceC*(cfb.dr(j-1)*(i-1))*(cfb.dr(j-1)*(i-1)));
                    end
                    

                    % dr2 = (cfb.dr(j-1)*(i-1))^2;
                    % quad  = quad + ((sourceC*dr2*cfb.dr(j-1))/sqrt(1-dr2));

                end

                quad = (cfb.dr(j-1)/2)*(left+right);
                cfb.R(j) = ((3*cfb.dt(j-1)*cfb.mu(cfb.t(j-1))) / (cfb.R(j-1)^2)) * quad + cfb.R(j-1);
                cfb.dr(j) = cfb.R(j)/(cfb.rSteps-1); % check indexes of all of these this is probably minus 1
                cfb.dt(j) = cfb.dt(j-1);
                cfb.t(j) = cfb.t(j-1) + cfb.dt(j); 

                cfb.C'
                cfb.T'
                cfb.W'
            end

            t = cfb.t;
            R = cfb.R;

        end

        function c = C_IC(cfb, r)
            d = 0.1; % parameterize later
            c = 0.9*cfb.K*((max(0, r - cfb.R(1) + d)^2)/(d^2));
        end
        
        function t = T_IC(cfb, r)
            d = 0.1; % parameterize later
            t = cfb.T0*min(((r - cfb.R(1) + d)^2)/(d^2), 1);
        end
        
        function w = W_IC(cfb, r)
            d = 0.1; % parameterize later
            w = (cfb.W0-cfb.Wh)*((max(0, r - cfb.R(1) + d)^2)/(d^2)); % look into w0
        end
        
        function c = C_BC(~, C, alpha, dr)
            c = (dr-(dr/alpha)+1)*C;
        end
        
        function t = T_BC(cfb, T, T0, dr)
            t = -cfb.gamma(T, T0)*dr*(T-T0) + T;
        end
        
        function w = W_BC(cfb, t)
            w = cfb.W0; % W ? what is this
        end

        function m = mu(~, t)
            m = 2 + (5/t);
        end

        function g = gamma(~, T, T0)
            if T <= T0 
                g = 1;
            else
                g = 0;
            end
        end
    end
end
