initF1 = mean_ft;
[LL,Res]=kf_tvp_corr1(FullParams,y, initF1);

% initF1 = mean_ft(1:3);
% [LL,Res]=kf_tvp_variance_only(EstimParams,y, initF1);

[a_sm, v_sm] = kf_smooth_original(Res.v_t, Res.invF_t, Res.L_t, Res.alfa_t(:,1:(end-1)), Res.P_t(:,:,1:(end-1)), Res.Z);
color = "g";
alpha=0.3;
fig = figure('Visible', 'off');
% set(fig, 'Position', [100, 100, 1600, 4]);
plot(1:size(y,2), y, 'LineWidth', 1);
hold on
confplot(1:(size(y,2)-2), a_sm(1,2:end),squeeze(sqrt(v_sm(1,1,2:end)) * 1.96)', color, alpha);
legend('px', 'nav','state');
hold off
saveas(fig, strcat('px_nav_state_plot/',"trial1",'.png'));



fig = figure('Visible', 'off');
% set(fig, 'Position', [100, 100, 1600, 4]);
plot(1:(size(y,2)-1), y(:,2:end), 'LineWidth', 1);
hold on
confplot(1:(size(y,2)-1), Res.alfa_t(1,2:end),squeeze(sqrt(Res.P_t(1,1,2:end)) * 1.96)', color, alpha);
legend('px', 'nav','state');
hold off
saveas(fig, strcat('px_nav_state_plot/',"trial1",'.png'));


% state_res = [a_sm(1,:);squeeze(sqrt(v_sm(1,1,:)))'];
% save(strcat('../state_result/',px_nav.isin{j}), 'state_res');