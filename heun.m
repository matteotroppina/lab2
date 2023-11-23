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
len =length(y_0);
y = zeros(len,t_size);
y(:,1) = y_0;

% check whether the input function is vector valued function 
if len > 1
    % explicit Heun
    for i = 1:t_size-1
        k1 = f(t_vector(i), y(:,i));
        k2 = f(t_vector(i+1), y(:,i) + dt * k1);
        y(:,i+1) = y(:,i) + dt/2 * (k1 + k2);
    end
    % Plot each component of f separately X vs t and Y vs t
    subplot(2,2,1);               
    plot(t_vector, y(1,:), 'b',"LineWidth",1.5)
    xlabel('t');
    ylabel('x(t)');
    subplot(2,2,2);               
    plot(t_vector, y(2,:), 'r',"LineWidth",1.5)
    xlabel('t');
    ylabel('y(t)');
    % plot y vs x in Van der Pol 
    subplot(2,2,[3,4]);
    plot(y(1,:) , y(2,:), 'k',"LineWidth",1.5)
    xlabel('x');
    ylabel('y');
    grid on;
else 
    for i = 1:t_size-1
        k1 = f(t_vector(i), y(i));
        k2 = f(t_vector(i+1), y(i) + dt * k1);
        y(i+1) = y(i) + dt/2 * (k1 + k2);
    end
    plot(t_vector, y, '-*'); 
end

end
