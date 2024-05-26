run Useful_Transformations;
px_nav_tot = readtable("mydata/px_nav.csv");
isin_list = unique(px_nav_tot.isin);

subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

for i = 1:length(isin_list)
    subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
end

for i=1:length(isin_list)
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
    
    px_nav.ldp = [NaN; px_nav.dprice(1:end-1)];
    px_nav.ldn = [NaN; px_nav.dnav(1:end-1)];
    px_nav.ptilde = px_nav.dprice-psi*px_nav.ldp;
    px_nav.ntilde = px_nav.dnav-phi*px_nav.ldn;
    
    knn = nanmean(px_nav.ntilde.^2);
    knp = nanmean(px_nav.ntilde.*px_nav.ptilde);
    kpp = nanmean(px_nav.ptilde.^2);
    knln = nanmean(px_nav.ntilde.*[NaN;px_nav.ntilde(1:end-1)]);
    knlp = nanmean(px_nav.ntilde.*[NaN;px_nav.ptilde(1:end-1)]);
    kpln = nanmean(px_nav.ptilde.*[NaN;px_nav.ntilde(1:end-1)]);
    kplp = nanmean(px_nav.ptilde.*[NaN;px_nav.ptilde(1:end-1)]);
    
    
    objfunc = @(vparam) moment_rec(vparam,[knn, kpp, knp, knln,...
    kplp, knlp, kpln], [psi, phi]);
    x0 = [0,0,0,0,0,0];
    lb = [1e-10, 1e-10, 1e-10, -Inf, -Inf, -Inf];
    ub = [Inf, Inf, Inf, Inf, Inf, Inf];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = @nonlcon;
    problem = createOptimProblem('fmincon', ...
    'objective', objfunc, ...
    'x0', x0, ...
    'lb', lb, ...
    'ub', ub, ...
    'nonlcon', @nonlcon,...
    'options', optimoptions('fmincon', 'Algorithm', 'sqp','TolFun',1e-12, 'TolX', 1e-7, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000, 'TolCon', 1e-20));
    gs = GlobalSearch;
    [t] = evalc('[x, fval, flag, out, allmins] = run(gs, problem);');
    
    cov_init = [x(1) x(4) x(5); x(4) x(2) x(6); x(5) x(6) x(3)];
    HL_init = chol(cov_init)';

    y = table2array(px_nav(:,{'lprice','lnav'}))';
    InitialParams = [InvUnoMenoUno(psi),InvUnoMenoUno(phi),HL_init(1,1), HL_init(2,2), HL_init(3,3), HL_init(2,1), HL_init(3,1), HL_init(3,2)];
    LossToMinmize = @(vparam) -base_kf_missing(vparam,y)/(size(y,2)-1);
    optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
    [EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);
    
    names = {"psi", "phi", "H11", "H22", "H33", "H21", "H31", "H32"};
    param_tab = table(names(:), EstimParams', 'VariableNames', {'Names', 'Value'});
    save(strcat('param_est/',px_nav.isin{1}), 'param_tab');
    
    % fit with optimal params and get smoothed states
    [LL,Res] = base_kf_missing(EstimParams,y);
    [a_sm, v_sm] = kf_smooth(Res.v_t, Res.invF_t, Res.L_t, Res.alfa_t(:,1:(end-1)), Res.P_t(:,:,1:(end-1)), Res.Z);
    color = [0.5, 0.5, 0.5];
    alpha=0.3;
    fig = figure;
    set(fig, 'Position', [100, 100, 1200, 800]);
    plot(px_nav.date, y);
    hold on
    confplot(px_nav.date, a_sm(1,:),squeeze(sqrt(v_sm(1,1,:)) * 1.96)', color, alpha);
    legend('px', 'nav','state');
    hold off
    saveas(fig, strcat('px_nav_state_plot/',px_nav.isin{1},'.png'));
    
    state_res = [a_sm(1,:);squeeze(sqrt(v_sm(1,1,:)))'];
    save(strcat('state_result/',px_nav.isin{1}), 'state_res');
    
    residual = cell(length(Res.v_t));
    for i=1:length(Res.v_t)
    residual{i} = chol(Res.invF_t{i})*Res.v_t{i};
    end
    
    
    for tt=1:length(residual)
    yt = y(:,tt);
    [Wt,nt]=SelectMatW(yt);
    if Wt==[1 0]
    residual{tt} = [residual{tt}; NaN];
    elseif Wt==[0 1]
    residual{tt} = [NaN; residual{tt}];
    elseif Wt==eye(2)
    residual{tt} = residual{tt};
    else
    residual{tt} = [NaN; NaN];
    end
    end

    % Diagnostics Check
    resid = [residual{:}];
    gen_diagnostics(px_nav, resid)
end
