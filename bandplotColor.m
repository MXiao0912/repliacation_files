function bandplotColor(x, y_bounds, color, alpha)
    % Ensure x and y_bounds are column vectors
    x = x(:);
    y_bounds = y_bounds';

    % Prepare data for the shaded area
    x_combined = [x; flipud(x)];
    y_combined = [y_bounds(1, :), flipud(y_bounds(2, :))];

    % Plot the shaded area
    fill(x_combined, y_combined, color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end