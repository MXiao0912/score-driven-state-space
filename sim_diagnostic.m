% get residual for diagnostic tests
residual = cell(length(Res.v_t));
for i=1:length(Res.v_t)
    residual{i} = chol(Res.invF_t{i})*Res.v_t{i};
end

for tt=1:length(residual)
    yt = y(:,tt+1)-y(:,tt);
    [Wt,nt]=SelectMatW(yt);
    if isequal(Wt,[1 0])
        residual{tt} = [residual{tt}; NaN];
    elseif isequal(Wt,[0 1])
        residual{tt} = [NaN; residual{tt}];
    elseif isequal(Wt,eye(2))
        residual{tt} = residual{tt};
    else
        residual{tt} = [NaN; NaN];
    end
end

% Diagnostics Check
resid = [residual{:}];
[p_jb, p_arch, p_lb]=gen_diagnostics(obs, resid);
% names = {"isin","p_jb", "p_arch", "p_lb"}; 
% test_tab = table(names(:), {"trail", p_jb, p_arch, p_lb}', 'VariableNames', {'Names', 'Value'});
% save(strcat('test/',"trail"), 'test_tab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot estimated ft
fig = figure('Visible', 'off');
subplot(4,1,1);
plot(px_nav.date, squeeze(Res.Q(1,1,2:end)));
hold on;
plot(px_nav.date, exp(ft(1,3:end)).^2);
subplot(4,1,2);
plot(px_nav.date, squeeze(Res.Q(2,2,2:end)));
hold on;
plot(px_nav.date, exp(ft(2,3:end)).^2);
subplot(4,1,3);
plot(px_nav.date, squeeze(Res.Q(3,3,2:end)));
hold on;
plot(px_nav.date, exp(ft(3,3:end)).^2);
subplot(4,1,4);
plot(px_nav.date, squeeze(Res.Q(2,1,2:end)./(Res.Q(1,1,2:end).*Res.Q(2,2,2:end)).^(0.5)));
hold on;
plot(px_nav.date, tanh(ft(4,3:end)));
saveas(fig, strcat('Diagnostics/trail_ft','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot score
% Res.score(:,500) = zeros(3,1)
fig = figure('Visible', 'off');
subplot(6,1,1);
plot(px_nav.date, squeeze(Res.score(1,:)));
subplot(6,1,2);
plot(px_nav.date, squeeze(Res.score(2,:)));
subplot(6,1,3);
plot(px_nav.date, squeeze(Res.score(3,:)));
subplot(6,1,4);
plot(px_nav.date, squeeze(Res.score(4,:)));
subplot(6,1,5);
plot(px_nav.date, squeeze(Res.score(5,:)));
subplot(6,1,6);
plot(px_nav.date, squeeze(Res.score(6,:)));
saveas(fig, strcat('trail_score','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot grad
fig = figure('Visible', 'off');
subplot(5,1,1);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.GradD{i}(1), 1:size(Res.GradD,1)));
subplot(5,1,2);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.GradD{i}(2), 1:size(Res.GradD,1)))
subplot(5,1,3);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.GradD{i}(3), 1:size(Res.GradD,1)))
% subplot(5,1,4);
% plot(1:size(Res.GradD,1), arrayfun(@(i) Res.GradD{i}(4), 1:size(Res.GradD,1)))
% subplot(5,1,5);
% plot(1:size(Res.GradD,1), arrayfun(@(i) Res.GradD{i}(5), 1:size(Res.GradD,1)))
saveas(fig, strcat('GRAD','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot info elements
fig = figure('Visible', 'off');
subplot(3,1,1);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.InfMat_temp{i}(1,1), 1:size(Res.GradD,1)))
subplot(3,1,2);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.InfMat_temp{i}(2,2), 1:size(Res.GradD,1)))
subplot(3,1,3);
plot(1:size(Res.GradD,1), arrayfun(@(i) Res.InfMat_temp{i}(3,3), 1:size(Res.GradD,1)))
saveas(fig, strcat('info_ele','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot infmat
fig = figure('Visible', 'off');
plot(1:size(Res.InfMat_temp,1), cellfun(@det, Res.InfMat_temp))
saveas(fig, strcat('infmat_temp','.png'));

fig = figure('Visible', 'off');
plot(1:size(Res.InfMat,1), cellfun(@det, Res.InfMat))
saveas(fig, strcat('infmat','.png'));

% plot invF
fig = figure();
plot(1:size(Res.invF_t,1), arrayfun(@(i) Res.invF_t{i}(1,1), 1:size(Res.invF_t,1)))
saveas(fig, "temp.png");

% plot vt
fig = figure();
plot(1:size(Res.v_t,1), arrayfun(@(i) Res.v_t{i}(1), 1:size(Res.v_t,1)))
saveas(fig, "temp.png");

% plot standardized vt
fig = figure();
plot(1:size(Res.v_t,1), arrayfun(@(i) sqrt(Res.invF_t{i}(1,1))*Res.v_t{i}(1), 1:size(Res.v_t,1)))
saveas(fig, "temp.png");