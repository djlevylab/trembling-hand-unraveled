function [PredictorsForRegression, PredictorsNamesForRegression, winningModelLambda, winningModelAlpha, winningModelMSE] = elastic_net(mouse_features,MMI,titles_for_elastic_net)

% Elastic net
alpha_val=0.1:0.1:1;

%run elastic net for selecting most important features
for iteration = 1:length(alpha_val)
    
    alpha=alpha_val(iteration);
    [B, FitInfo] = lasso(mouse_features,MMI,'CV',10,'alpha',alpha,'PredictorNames',titles_for_elastic_net,'Standardize',false);
    idxLambdaMinMSE(iteration,1) = FitInfo.IndexMinMSE(1);
    MinMSE(iteration,1) = FitInfo.MSE(idxLambdaMinMSE(1));
    Lambda(iteration,1) = FitInfo.LambdaMinMSE(1);
    minMSEModelPredictorNames{iteration,1} = FitInfo.PredictorNames(B(:,idxLambdaMinMSE(1))~=0);
  
    clear B FitInfo
    
end

% Choose the winning elastic net model (min MSE), report alpha and lambda

[winningModelMSE, idx_winning_model] = min(MinMSE);
winningModelLambda = Lambda(idx_winning_model);
winningModelAlpha = alpha_val(idx_winning_model);
winningModelPredictorNames = minMSEModelPredictorNames{idx_winning_model,1};
winningModelPredictors_logical = ismember(titles_for_elastic_net,winningModelPredictorNames);
winningModelPredictors_idx = find(winningModelPredictors_logical==1);

PredictorsForRegression = mouse_features(:,winningModelPredictors_idx);
PredictorsNamesForRegression = titles_for_elastic_net(winningModelPredictors_idx);

end