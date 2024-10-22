EstimParams_all = NaN(14,100);
for i=1:100
    EstimParams = load("verify_var_res/EstimParams_res_" + i + ".mat");
    EstimParams = EstimParams.EstimParams;
    EstimParams_all(:,i) = EstimParams;
end