function[] = gen_diagnostics(px_nav, data)
% residual plot
fig = figure;
set(fig, 'Position', [100, 100, 1200, 800]);

subplot(3,2,1)
plot(px_nav.date, data(1,:));
title("standardized residuals - lprice");
subplot(3,2,2)
plot(px_nav.date, data(2,:));
title("standardized residuals - lnav");
% autocorrelation plots
subplot(3,2,3);
plot_acf(data(1,:));
% Add title and labels
title('ACF lprice');
xlabel('Lag');
ylabel('Autocorrelation');
subplot(3,2,4);
plot_acf(data(2,:));
% Add title and labels
title('ACF lnav');
xlabel('Lag');
ylabel('Autocorrelation');
% tests
% Normality Test (Jarque-Bera)
[~, p_jb, jb_stat] = jbtest(data(1,:));
% Heteroskedasticity Test (ARCH)
[~, p_arch, arch_stat] = archtest(fillmissing(data(1,:), "linear"));
% Autocorrelation Test (Ljung-Box Q-test) 20 lags
[~, p_lb, lb_stat] = lbqtest(data(1,:));
% Create a table to summarize the results
testNames = {'Jarque-Bera Test'; 'ARCH Test'; 'Ljung-Box Q-test'};
pValues = [p_jb; p_arch; p_lb];
testStats = [jb_stat; arch_stat; lb_stat];
resultsTable = table(testNames, pValues, testStats, ...
'VariableNames', {'Test', 'PValue', 'TestStatistic'});
% Plot 3: Test Results Text
subplot(3, 2, 5);
axis off; % Turn off the axis
% Prepare the text for display
resultsText = {
sprintf('Test: %s', 'Jarque-Bera Test');
sprintf('P-Value: %.4f', p_jb);
sprintf('Test Statistic: %.4f', jb_stat);
'';
sprintf('Test: %s', 'ARCH Test');
sprintf('P-Value: %.4f', p_arch);
sprintf('Test Statistic: %.4f', arch_stat);
'';
sprintf('Test: %s', 'Ljung-Box Q-test');
sprintf('P-Value: %.4f', p_lb);
sprintf('Test Statistic: %.4f', lb_stat);
};
% Display the text in the subplot
text(0, 0.5, resultsText, 'FontName', 'Courier', 'FontSize', 10);
title('Test Results -- lprice');
% Normality Test (Jarque-Bera)
[~, p_jb, jb_stat] = jbtest(data(2,:));
% Heteroskedasticity Test (ARCH)
[~, p_arch, arch_stat] = archtest(fillmissing(data(2,:), "linear"));
% Autocorrelation Test (Ljung-Box Q-test) 20 lags
[~, p_lb, lb_stat] = lbqtest(data(2,:));
% Create a table to summarize the results
testNames = {'Jarque-Bera Test'; 'ARCH Test'; 'Ljung-Box Q-test'};
pValues = [p_jb; p_arch; p_lb];
testStats = [jb_stat; arch_stat; lb_stat];
resultsTable = table(testNames, pValues, testStats, ...
'VariableNames', {'Test', 'PValue', 'TestStatistic'});
% Plot 3: Test Results Text
subplot(3, 2, 6);
axis off; % Turn off the axis
% Prepare the text for display
resultsText = {
sprintf('Test: %s', 'Jarque-Bera Test');
sprintf('P-Value: %.4f', p_jb);
sprintf('Test Statistic: %.4f', jb_stat);
'';
sprintf('Test: %s', 'ARCH Test');
sprintf('P-Value: %.4f', p_arch);
sprintf('Test Statistic: %.4f', arch_stat);
'';
sprintf('Test: %s', 'Ljung-Box Q-test');
sprintf('P-Value: %.4f', p_lb);
sprintf('Test Statistic: %.4f', lb_stat);
};
% Display the text in the subplot
text(0, 0.5, resultsText, 'FontName', 'Courier', 'FontSize', 10);
title('Test Results -- lnav');
% Adjust the layout for better visualization
sgtitle('Diagnostics Checks');

saveas(fig, strcat('Diagnostics_base_individual/',px_nav.isin{1},'Diagnostics.png'));
end