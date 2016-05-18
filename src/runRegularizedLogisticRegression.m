function [ results ] = runRegularizedLogisticRegression( X,y, testIdx, options )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% hold out the test set
X_train = X(~testIdx,:);
X_test = X(testIdx,:);
y_train = y(~testIdx);
y_test = y(testIdx);

% fit lasso
options.alpha = 1; % 1 == lasso, 0 == ridge
cvfit = cvglmnet(X_train, y_train, 'binomial', options);
% save coeff
% results.lasso.coef_1se = cvglmnetCoef(cvfit, 'lambda_1se');
% results.lasso.lambda_1se = cvfit.lambda_1se;
% results.lasso.coef_min = cvglmnetCoef(cvfit, 'lambda_min');
% results.lasso.lambda_min = cvfit.lambda_min;

% compute the performance
y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_1se));
results.lasso.accuracy.onese = sum(y_hat == y_test) / length(y_test);
y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_min));
results.lasso.accuracy.min = sum(y_hat == y_test) / length(y_test);

% % fit ridge
% options.alpha = 0; % 1 == lasso, 0 == ridge
% cvfit = cvglmnet(X_train, y_train, 'binomial', options);
% % save coeff
% results.ridge.coef_1se = cvglmnetCoef(cvfit, 'lambda_1se');
% results.ridge.lambda_1se = cvfit.lambda_1se;
% results.ridge.coef_min = cvglmnetCoef(cvfit, 'lambda_min');
% results.ridge.lambda_min = cvfit.lambda_min;
% 
% % compute the performance
% y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_1se));
% results.ridge.accuracy.onese = sum(y_hat == y_test) / length(y_test);
% y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_min));
% results.ridge.accuracy.min = sum(y_hat == y_test) / length(y_test);

end

