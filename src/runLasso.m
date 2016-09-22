function [ results ] = runLasso(X, y, testIdx, options)
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
results.lasso_lambda_min = cvfit.lambda_min;
results.lasso_coef_lambda_min = cvglmnetCoef(cvfit, 'lambda_min');

% compute the performance
y_hat = myStepFunction(cvglmnetPredict(cvfit, X_test,cvfit.lambda_min));
results.lasso_accuracy_lambda_min = sum(y_hat == y_test) / length(y_test);
end

