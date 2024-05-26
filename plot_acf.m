function[] = plot_acf(data)
    % Define the number of lags
    numLags = 20; % Number of lags
    
    % Calculate the ACF values
    [acfValues, lags, bounds] = autocorr(data, 'NumLags', numLags);
    
    % Exclude the 0 lag (first element)
    acfValues = acfValues(2:end);
    lags = lags(2:end);
    
    % Plot the ACF values excluding 0 lag
    stem(lags, acfValues, 'filled');
    hold on;
    
   % Define the x coordinates for the polygon
    x = [1, numLags, numLags, 1];
    
    % Ensure y coordinates are properly set
    y = [-bounds(1), -bounds(1), bounds(1), bounds(1)];

    % Plot the shaded area for the confidence bounds
    fill(x, y, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    hold off;

end
