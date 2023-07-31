function [data] = logistic_regression_data_generator(n, d, alpha)
% Data generator for logistic regression problem.
%
% Inputs:
%       n               number of samples.
%       d               number of dimensions.
%       alpha           1 x n vector specifying Lipschitz weights
% Output:
%       data            data set
%       data.x_train    train data of x.
%       data.y_train    train data of y.
%       data.x_test     test data of x.
%       data.y_test     test data of y.
%       data.w_opt      solusion.
%
% This file is part of SGDLibrary.
%
% Created H.Kasai on Oct. 25, 2016

% PLEASE NOTE: THIS IS A SINGLE FILE FROM ANOTHER THE SGDLibrary.
% THE ORIGINAL REPO WITH THE FULL LIBRARY (NOT OWNED BY ME!) IS AVAILABLE
% AT
% https://github.com/hiroyuki-kasai/SGDLibrary
% AND SEE THE CORRESPONDING PAPER
% Kasai, Hiroyuki. "SGDLibrary: A MATLAB library for stochastic 
% optimization algorithms." 
% The Journal of Machine Learning Research 18.1 (2017): 7942-7946.

    w_opt = randn(d, 1);  
    %w_opt = 0.5 * ones(d, 1);
    data.w_opt = w_opt;    

    % train data
    x1 = repmat(alpha,d,1) .* randn(d, n);
    y1 = rand(1, n) < sigmoid(w_opt' * x1);
    y1 = 2*y1 - 1;
    assert(sum(y1 == 1) + sum(y1 == -1) == n);

    data.x_train = x1;
    data.y_train = y1;
    
    % test data    
    x2 = 20 * randn(d, n);
    y2 = rand(1, n) < sigmoid(w_opt' * x2);
    y2 = 2*y2 - 1;
    assert(sum(y2 == 1) + sum(y2 == -1) == n);
    
    data.x_test = x2;
    data.y_test = y2;

    data.w_init = w_opt + (1/max(alpha))*randn(d,1);
end

function [y] = sigmoid(x)
% Sigmoid function.

    y = 1.0 ./ (1.0 + exp(-x));

end