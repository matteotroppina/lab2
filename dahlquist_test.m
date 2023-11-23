function error_matrix = dahlquist_test(lambda, x_0, t_end, dt_vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dahlquistz_test - y(t) = lambda * x
%
%   dahlquist_test(x_0, t_end, dt_vector, f, f_sol) performs simulations of
%   various numerical methods (Euler, Heun, RK4) for solving ODEs in the form
%
%                       y(t) = lambda * x
%
%   and analyzes their errors.
%
%   Input:
%       lambda      - y(t) = lambda * x
%       x_0         - initial value
%       t_end       - final time
%       dt_vector   - time step sizes vector

%   Output:
%       error_matrix    - 3 x N matrix (N = number of time step sizes)
%                         each column corresponds to a time step size while
%                         the rows to the three methods implemented
%
%   This function generates multiple plots to visualize the results of the
%   simulations, including the analytical solution, Euler method, Heun method,
%   and RK4 method. It also calculates and displays errors for each time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = @(t,x) lambda*x; %ODE
f_sol = @(t) exp(lambda*t); %ODE analytical solution

%vector for x axis plot
x_vector = linspace(0,t_end,501);

%error function
error_matrix = zeros(3,length(dt_vector));
E = @(y, dt, t_end, x_exact) sqrt( (dt/t_end) * sum( (y - x_exact).^2 ));

%legend strings
legend_str = strings(1, length(dt_vector) + 1);
legend_str(1) = 'y(t) = \ite^{\lambda t}';


% analytical solution plot
figure(1)
plot(x_vector, f_sol(x_vector), 'k', "LineWidth", 1.5);
title('Analytical Solution', legend_str(1));
xlabel('t')
ylabel('y(t)')

% Euler method figure
figure(2);
plot(x_vector, f_sol(x_vector), 'k', "LineWidth", 1.5);
title('Euler Method');
xlabel('t')
ylabel('y(t)')

%Heun method figure
figure(3);
plot(x_vector, f_sol(x_vector), 'k', "LineWidth", 1.5);
title('Heun Method');
xlabel('t')
ylabel('y(t)')

%rk4 method figure
figure(4);
plot(x_vector, f_sol(x_vector), 'k', "LineWidth", 1.5);
title('RK4 Method');
xlabel('t')
ylabel('y(t)')

for i = 1:length(dt_vector)

    dt = dt_vector(i);
    t_vector = 0:dt:t_end;    %time steps vector
    legend_str(i + 1) = sprintf('\\deltat_{%d} = %.3f', i, dt); %fills the legend strings

    %expl_euler
    figure(2);
    hold on;
    y_euler = expl_euler(x_0, dt, t_end, f);
    error_matrix(1,i) = E(y_euler, dt, t_end, f_sol(t_vector)); %errors for euler are stored in the first row
    text(i, 0.55, sprintf('e_{%d} = \n %.4e', i, error_matrix(1,i))); %prints errors on plot

    %heun
    figure(3);
    hold on;
    y_heun = heun(x_0, dt, t_end, f);
    error_matrix(2,i) = E(y_heun, dt, t_end, f_sol(t_vector)); %errors for heun are stored in the second row
    text(i, 0.55, sprintf('e_{%d} = \n %.4e', i, error_matrix(2,i))); %prints errors on plot
    
    %rk4
    figure(4);
    hold on;
    y_rk4 = runge_kutta_4(x_0, dt, t_end, f);
    error_matrix(3,i) = E(y_rk4, dt, t_end, f_sol(t_vector)); %errors for rk4 are stored in the third row
    text(i, 0.55, sprintf('e_{%d} = \n %.4e', i, error_matrix(3,i))); %prints errors on plot

end

figure(1)
legend(legend_str(1));

figure(2)
legend(legend_str);

figure(3)
legend(legend_str);

figure(4)
legend(legend_str);

end
