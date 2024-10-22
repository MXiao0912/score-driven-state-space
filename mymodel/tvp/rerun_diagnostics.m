function run_pilot_mymodel(j)
    j = str2double(j);
    % j=12;
    cd("/rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/")
    addpath("../")
    run Useful_Transformations.m

    px_nav_tot = readtable("../../../mydata/px_nav.csv");
    isin_list = unique(px_nav_tot.isin);

    subTables = cell(length(isin_list), 1);  % Preallocate a cell array to hold sub-tables

    for i = 1:length(isin_list)
        subTables{i} = px_nav_tot(strcmp(px_nav_tot.isin,isin_list(i)), :);
    end

    px_nav = subTables{j};
    px_nav.lprice = log(px_nav.price);
    px_nav.lnav = log(px_nav.nav);
    y = table2array(px_nav(:,{'lprice','lnav'}))';
    
    Res = load(strcat('Resall/',px_nav.isin{j})).Res;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic checks
    residual = cell(length(Res.v_t));
    for i=1:length(Res.v_t)
        try
            residual{i} = chol(Res.invF_t{i})*Res.v_t{i};
        catch
            % Code to execute if the try block fails
            disp('An error occurred, executing alternative code...');
            [Wt,nt]=SelectMatW(y(:,i));
            if isequal(Wt,[1 0])
                residual{tt} = NaN;
            elseif isequal(Wt,[0 1])
                residual{tt} = NaN;
            elseif isequal(Wt,eye(2))
                residual{tt} = [NaN; NaN];
            else
                residual{tt} = [NaN; NaN];
            end
        end
    end

    for tt=1:length(residual)
        [Wt,nt]=SelectMatW(y(:,tt));
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
    [p_jb_p, p_arch_p, p_lb_p,p_jb_n, p_arch_n, p_lb_n]=gen_diagnostics(px_nav, resid);
    names = {"isin","p_jb_p", "p_arch_p", "p_lb_p","p_jb_n", "p_arch_n", "p_lb_n"}; 
    test_tab = table(names(:), {isin_list{j}, p_jb_p, p_arch_p, p_lb_p,p_jb_n, p_arch_n, p_lb_n}', 'VariableNames', {'Names', 'Value'});
    save(strcat('test/',px_nav.isin{j}), 'test_tab');
end