%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating and visualizing Chebyshev points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a function handle for the Chebyshev Polynomial of
% 2nd kind, and its roots (which are Chebyshev points)
U       = @(n, x) cos(n * acos(x));
P       = @(n) cos([n:-1:-n] * pi / n);
order   = 0;

% Sample points at which we will evaluate Tn(x)
npts    = 1000;
pts     = sort(2*rand(npts,1) - 1);

figure;
orders_to_try = 2.^[3:4]';
n_orders      = length(orders_to_try)
for order = 1:n_orders
    vals  = U(orders_to_try(order), pts);
    subplot(n_orders, 1, order);
    cheb_pts = P(orders_to_try(order));
    scatter(cheb_pts, zeros(size(cheb_pts)), 'MarkerFaceColor', 'b');
    hold on;
    plot(pts, vals, '-.', 'LineWidth', 2.0);
    xlabel('x');
    ylabel(['T_{', int2str(orders_to_try(order)), '}(x)'])
    set(gca, 'FontSize', 28)
    grid on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
