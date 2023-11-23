function y = heun(y_0, dt, t_end, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heun - implements the Heun method for solving ordinary
%        differential equations (ODEs).
%
%        Heun method:
%        y_n+1 = y_n + Delta_t * 1/2*( dy_n(t_n,y_n) + 
%                                      dy_n+1(t_n+1, y_n + dt*dy_n(t_n,y_n) ) 
%       
% y = heun(y_0, dt, t_end, f)
%
% Inputs:
%       y_0   - initial condition
%       dt    - time step size for the integration
%       t_end - final time
%       f     - function handle of the ODE, must be of the form:
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

% explicit Heun
for i = 1:t_size-1
    y(i+1) = y(i) + dt/2 * ( f(t_vector(i), y(i)) + ...
                             f(t_vector(i+1), y(i) + dt * f(t_vector(i), y(i)) ) );
end

plot(t_vector, y, '-*');    %plot y(t)


end
