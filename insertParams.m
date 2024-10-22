function params = insertParams(variableParams, fixedParams)
    params = fixedParams;
    varIdx = 1;
    for i = 1:length(params)
        if isnan(params(i))
            params(i) = variableParams(varIdx);
            varIdx = varIdx + 1;
        end
    end
end