cd("../baseline_original/code/")
[psi, phi, x] = nonparam(obs);
    
% kalman filter estimation
y = table2array(obs(:,{'price','nav'}))';
InitialParams = [InvUnoMenoUno(psi),InvUnoMenoUno(phi),0,0,0,sqrt(x)];
LossToMinmize = @(vparam) -base_kf_original(vparam,y)/(size(y,2)-1);
optionsIVAN = optimset('Display', 'iter-detailed','LargeScale', 'off','MaxFunEvals',5000);
[EstimParams] = fminunc(LossToMinmize,InitialParams,optionsIVAN);

names = {"psi", "phi", "mean_p", "mean_n", "mean_r", "H11", "H22", "H33"}; % H here is unrecovered variances +-sqrt(.)
param_tab = table(names(:), EstimParams', 'VariableNames', {'Names', 'Value'});
save(strcat('../param_est/',"trail"), 'param_tab');

% fit with optimal params and get smoothed states
[LL,Res] = base_kf_original(EstimParams,y);
[a_sm, v_sm] = kf_smooth_original(Res.v_t, Res.invF_t, Res.L_t, Res.alfa_t(:,1:(end-1)), Res.P_t(:,:,1:(end-1)), Res.Z);
color = "g";
alpha=0.3;
fig = figure('Visible', 'off');
set(fig, 'Position', [100, 100, 1600, 400]);
plot(1:size(y,2), y, 'LineWidth', 1);
hold on
confplot(1:size(a_sm(1,2:end),2), a_sm(1,2:end),squeeze(sqrt(v_sm(1,1,2:end)) * 1.96)', color, alpha);
legend('px', 'nav','state');
hold off
saveas(fig, strcat('../px_nav_state_plot/',"trail",'.png'));

% state_res = [a_sm(1,:);squeeze(sqrt(v_sm(1,1,:)))'];
% save(strcat('../state_result/',px_nav.isin{j}), 'state_res');

% % get residual for diagnostic tests
% residual = cell(length(Res.v_t));
% for i=1:length(Res.v_t)
%     residual{i} = chol(Res.invF_t{i})*Res.v_t{i};
% end

% for tt=1:length(residual)
% yt = y(:,tt+1)-y(:,tt);
% [Wt,nt]=SelectMatW(yt);
% if isequal(Wt,[1 0])
% residual{tt} = [residual{tt}; NaN];
% elseif isequal(Wt,[0 1])
% residual{tt} = [NaN; residual{tt}];
% elseif isequal(Wt,eye(2))
% residual{tt} = residual{tt};
% else
% residual{tt} = [NaN; NaN];
% end
% end

% % Diagnostics Check
% resid = [residual{:}];
% [p_jb, p_arch, p_lb]=gen_diagnostics(px_nav(2:end,:), resid);
% names = {"isin","p_jb", "p_arch", "p_lb"}; 
% test_tab = table(names(:), {isin_list{j}, p_jb, p_arch, p_lb}', 'VariableNames', {'Names', 'Value'});
% save(strcat('test/',px_nav.isin{j}), 'test_tab');
