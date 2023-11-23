function van_oscillator(y_0, t_end, dt ,f)
figure(5);
hold on;
heun(y_0, dt, t_end, f);
hold off;

% create a mesh grid to calculate the derivative at each point in it
x = linspace(-5, 5, 20);
y = linspace(-5, 5, 20);
[xGrid, yGrid] = meshgrid(x, y);

% initialize values at each point in the grid
xGrid_values = zeros(size(xGrid));
yGrid_values = zeros(size(yGrid));

% calculate the derivative value at each point
for i = 1:numel(xGrid)
    dydt = f(0, [xGrid(i); yGrid(i)]);
    xGrid_values(i) = dydt(1);
    yGrid_values(i) = dydt(2);
end

figure(6);
quiver(xGrid, yGrid, xGrid_values, yGrid_values, 'AutoScale', 'on', 'LineWidth', 2);
xlabel('Displacement (x)');
ylabel('Velocity (dx/dt)');
title('Van der Pol Oscillator - Phase Space Quiver Plot');
grid on;

end