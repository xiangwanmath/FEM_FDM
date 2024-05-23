clc
clear

tSteps = 5;
rSteps = 100;
rInit = 0.5;

sim = FDM_1D(tSteps, rSteps, rInit);

[t, R] = sim.run_euler();

% disp(R);
% disp(t);

% plot(t, R);