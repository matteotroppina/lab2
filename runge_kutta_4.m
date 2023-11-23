function y = runge_kutta_4(y_0, dt, t_end, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runge_kutta_4 - implements the Runge-Kutta method of fourth order for solving 
%                 ordinary differential equations (ODEs)
%
%                 Runge Kutta method:
%                 y_n+1 = y_n + dt/6 * (Y_1 + 2 * Y_2 + 2 * Y_3 + Y_4);
%                 Y_1 = f (t_n, y_n);
%                 Y_2 = f (t_(n+1/2), y_n + dt/2 * Y_1);
%                 Y_3 = f (t_(n+1/2), y_n + dt/2 * Y_2);
%                 Y_4 = f (t_n+1, y_n + dt * Y_3);
%
% y = runge_kutta_4(y_0, dt, t_end, f)
%
% Inputs:
%       y_0   - initial condition
%       dt    - time step size for the integration
%       t_end - final time
%       f     - function handle representing the right-hand side of the ODE.
%               function must be of the form 
%
%                               dy = f(t, y)
%
%               where dy is the rate of change at time t, and y is the 
%               current value of the dependent variable.
%
% Outputs:
%       y     - vector of solution values. Each element represents the solution at
%               a specific time step, with the first row corresponding to y_0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vector of time steps (starts at t0)
t0 = 0;
t_vector = t0:dt:t_end;
t_size = length(t_vector);

%initialization of solutions vector
y = zeros(1,length(t_size));
y(1) = y_0;

% Runge Kutta (fourth order)
for i=1:t_size-1
    Y_1 = f (t_vector(i), y(i));
    Y_2 = f (t_vector(i) + dt/2, y(i) + dt/2 * Y_1);
    Y_3 = f (t_vector(i) + dt/2, y(i) + dt/2 * Y_2);
    Y_4 = f (t_vector(i) + dt, y(i) + dt * Y_3);

    y(i+1) = y(i) + dt/6 * (Y_1 + 2 * Y_2 + 2 * Y_3 + Y_4);
end

plot(t_vector, y, '-*');    %plot y(t)

end