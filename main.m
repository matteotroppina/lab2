clear
close all
clc

% lambda = -1;
% t_end = 5; % final time
% x_0 = 1; % intial condition
% dt_vector = [1, 0.5, 0.25, 0.125]; %vector of time steps size (dt)
% 
% f = @(t,x) lambda * x; %ODE
% f_sol = @(t) exp(lambda*t); %ODE analytical solution
% 
% dahlquist_test(x_0, t_end, dt_vector,f, f_sol);

% constants initialization
mu = 1;
t_end = 20;
dt = 0.1;
x_0 = 1;    % initial values of x
y_0 = 1;    % initial values of y

% declaring the function of Van der Pol Oscillator
f = @(t,y) [y(2) ; mu*(1-y(1)^2)*y(2) - y(1)];

van_oscillator([x_0,y_0], t_end, dt, f);

