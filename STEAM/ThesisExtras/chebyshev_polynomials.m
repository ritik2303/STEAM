%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating and visualizing Chebyshev polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a function handle for a Chebyshev Polynomial
T       = @(n, x) cos(n * acos(x));
order   = 0;

% Sample points at which we will evaluate Tn(x)
npts    = 1000;
pts     = sort(2*rand(npts,1) - 1);
vals    = T(order, pts);

% Use for a single plot
% figure;
% hold on;
% plot(pts, vals);
% xlabel('x');
% ylabel(['T_{',str(n), '}(x)'])

figure;
hold on; grid on;
orders_to_try = 2.^[0:3]';
n_orders      = length(orders_to_try);
lines         = zeros(n_orders, 1);
for order = 1:n_orders
    vals  = T(orders_to_try(order), pts);
    ln    = plot(pts, vals, 'LineWidth', 2.0);
    % text(pts(150), vals(150), ['n=',int2str(orders_to_try(order))], 'FontSize', 28);
    lines(order) = ln;
    xlabel('x');
    ylabel('T_{n}(x)')
end
base_string = repmat('n = ', size(orders_to_try));
legend(lines, [base_string, int2str(orders_to_try)]);
set(gca, 'FontSize', 28)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
