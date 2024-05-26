px_nav_tot = readtable("mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

for i=1:size(isin_list,1)
    px_nav = subTables{i};
    px_nav.lprice = log(px_nav.price);
    px_nav.lnav = log(px_nav.nav);
    px_nav.dprice = [NaN; diff(px_nav.lprice)];
    px_nav.dnav = [NaN; diff(px_nav.lnav)];
    px_nav.dprice = fillmissing(px_nav.dprice,'previous');
    px_nav.dnav = fillmissing(px_nav.dnav,'previous');
    
    model_arma_11 = arima('ARLags', 1, 'D', 0, 'MALags', 1);
    [estModelpx, estParamspx] = estimate(model_arma_11, px_nav.dprice);
    [estModelnav, estParamsnav] = estimate(model_arma_11, px_nav.dnav);
    psi = estModelpx.AR{1};
    phi = estModelnav.AR{1};
    t_ = (1/size(px_nav,1)):(1/size(px_nav,1)):1;
    silverman = @(n) 1.06*min(std(t_),iqr(t_)/1.34)/(n^0.2);
    h = silverman(size(px_nav,1));
    
    px_nav.ldp = [NaN; px_nav.dprice(1:end-1)];
    px_nav.ldn = [NaN; px_nav.dnav(1:end-1)];
    px_nav.ptilde = px_nav.dprice-psi*px_nav.ldp;
    px_nav.ntilde = px_nav.dnav-phi*px_nav.ldn;
    knn = smoothdata(px_nav.ntilde.^2, 'gaussian',h);
    knp = smoothdata(px_nav.ntilde.*px_nav.ptilde, 'gaussian',h);
    kpp = smoothdata(px_nav.ptilde.^2, 'gaussian',h);
    knln = smoothdata(px_nav.ntilde.*[NaN;px_nav.ntilde(1:end-1)], 'gaussian',h);
    knlp = smoothdata(px_nav.ntilde.*[NaN;px_nav.ptilde(1:end-1)], 'gaussian',h);
    kpln = smoothdata(px_nav.ptilde.*[NaN;px_nav.ntilde(1:end-1)], 'gaussian',h);
    kplp = smoothdata(px_nav.ptilde.*[NaN;px_nav.ptilde(1:end-1)], 'gaussian',h);

    figure;
    plot(px_nav.date,kplp)

    
    res = cell(1,size(knln,1));
    res_fval = cell(1,size(knln,1)); 
    res_flag = cell(1,size(knln,1));
    res_out = cell(1,size(knln,1));
    res_allmins = cell(1,size(knln,1));
    
    for j=4:size(knln,1)
        disp({i,j});
        objfunc = @(vparam) moment_rec(vparam,[knn(j,1), kpp(j,1), knp(j,1), knln(j,1),...
            kplp(j,1), knlp(j,1), kpln(j,1)], [psi, phi]);
        x0 = [0.01, 0.01, 0.01, 0];
        lb = [1e-10, 1e-10, 1e-10, -Inf];
        ub = [Inf, Inf, Inf, Inf];
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        nonlcon = @(x) deal(((x(4)^2)/(x(1)*x(3)))-1,[]);
        %[x, fval, exitflag, output]=fmincon(objfunc, x0, A, b, Aeq, beq, lb, ub, nonlcon, optimoptions('fmincon','Algorithm', 'sqp','TolFun',1e-20, 'TolCon',1e-20));
        %[x, fval, exitflag, output]=fmincon(objfunc, x0, A, b, Aeq, beq, lb, ub, nonlcon, optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-10));
        problem = createOptimProblem('fmincon', ...
            'objective', objfunc, ...
            'x0', x0, ...
            'lb', lb, ...
            'ub', ub, ...
            'nonlcon', nonlcon,...
            'options', optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-10));
        %ms = MultiStart('Display', 'iter', 'UseParallel', true, 'StartPointsToRun','bounds-ineqs');
        
        % Run MultiStart with fmincon
        %[x1, x2, x3, x4] = meshgrid(linspace(0.1,1,10), linspace(0.1,1,10), linspace(0.1,1,10), linspace(-0.9,0.9,20));
        %start_points = [x1(:),x2(:),x3(:),x4(:)];
        %start_points_set = CustomStartPointSet(start_points);
        %[x, fval, flag, out, allmins] = run(ms, problem, 100); % 50 runs from different start points
        gs = GlobalSearch;
        %[x, fval, flag, out, allmins] = run(gs, problem);
        [t] = evalc('[x, fval, flag, out, allmins] = run(gs, problem);');
    
        res{j}=x;
        res_fval{j}=fval;
        res_flag{j}=flag;
        res_out{j}=out;
        res_allmins{j}=allmins;
        
    end
    
    M = cell2mat(res);
    M = reshape(M,4,[]).';
    M(:,5) = M(:,4)./sqrt(M(:,1).*M(:,3));
    
    figure;  % Open a new figure window
    subplot(2,3,1);
    %plot(px_nav.date(2:end), abs(M(:,4)), 'o-', MarkerSize=1);  % Change 'o-' to any other plot style as needed
    plot(px_nav.date(4:end), M(:,5), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\rho_{\epsilon, \omega}');  % Label y-axis
    title('\rho_{\epsilon, \omega}');  % Add a title
    
    subplot(2,3,2);
    plot(px_nav.date(4:end), M(:,4), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('cov_{\epsilon, \omega}');  % Label y-axis
    title('cov_{\epsilon, \omega}');  % Add a title
    
    subplot(2,3,3);
    plot(px_nav.date(4:end), M(:,1), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_{\epsilon}^2');  % Label y-axis
    title('\sigma_{\epsilon}^2');  % Add a title
    
    subplot(2,3,4);
    plot(px_nav.date(4:end), M(:,2), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_r^2');  % Label y-axis
    title('\sigma_r^2');  % Add a title
    
    subplot(2,3,5);
    plot(px_nav.date(4:end), M(:,3), '-', MarkerSize=5);  % Change 'o-' to any other plot style as needed
    xlabel('dates');  % Label x-axis
    ylabel('\sigma_{\omega}^2');  % Label y-axis
    title('\sigma_{\omega}^2');  % Add a title
    
    sgtitle(strcat(unique(px_nav.isin),',', 'h=silverman'))
    
    saveas(gcf, strcat(px_nav.isin{1},',','h=silverman','.png'))


end 