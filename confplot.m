function confplot(x, y, e, color, alpha)
    % Ensure x, y, and e are column vectors
    x = x(:);
    y = y(:);
    e = e(:);

    % Calculate the upper and lower bounds
    upper_bound = y + e;
    lower_bound = y - e;

    % Combine the bounds into a matrix
    y_bounds = [lower_bound, upper_bound];

    % Ensure x and y_bounds are column vectors
    x = x(:);

    % Prepare data for the shaded area
    x_combined = [x; flipud(x)];
    y_combined = [y_bounds(:, 1); flipud(y_bounds(:, 2))];

    % Plot the shaded area using 'fill'
    fill(x_combined, y_combined, color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    
    % Optionally, plot the base of the shaded area (mm) if needed
    % This line is commented out to avoid showing the mm line
    % plot(x, repmat(mm, size(x)), 'Color', color, 'LineWidth', 1.5);

    % Overlay the main line
    hold on;
    plot(x, y, 'Color', color, 'LineWidth', 1.5);
    hold off;
end
