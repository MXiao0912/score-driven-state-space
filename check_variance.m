tot = NaN(14, 100);
for i = 1:100
    if i~=65
        EstimParamssim = load(sprintf("verify_var_res/EstimParams_res_%d.mat", i));
        tot(:,i) = EstimParamssim.EstimParams;
    end
end

tot_std = std(tot, 0, 2, 'omitnan');


tot_std =

  463.1669
    0.2357
    0.0064
    0.0162
    0.7074
    0.9575
    0.6828
    0.9692
    1.5899
    2.0429
    1.1952
    1.2638
    2.0598
    0.9858

>> check(2:end)./tot_std(2:end)

ans =

    0.0249
    1.1556
    0.2038
    0.0068
    0.3047
    0.2812
    1.3729
    0.4584
    3.3768
    5.6400
    2.0827
    4.5942
    0.3119
