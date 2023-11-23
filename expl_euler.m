function y = expl_euler(y_0, dt, t_end, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expl_euler - implements the explicit Euler method for solving ordinary
%              differential equations (ODEs). 
%
%              explicit Euler method:
%              % y_n+1 = y_n + Delta_t * dy_n(t_n,y_n)
%
% y = expl_euler(y_0, dt, t_end, f)
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

% explicit Euler
for i = 1:t_size-1
    y(i + 1) = y(i) + dt * f(t_vector(i), y(i));
end

plot(t_vector, y, '-*');    %plot y(t)


end


