function [F,J] =log_reg(x,Set,X,y,lambda)

% X is an n x m matrix, y is 1 x m vector

% F is 1 x n, J is m x n 
[n,m] = size(X);

%Set defines the component functions of F you would like 
if nargin<2 % Use all components
	Set = 1:length(x); %(in this case there are dim n components)
elseif length(Set) > length(unique(Set)) 
	disp('Warning: Set has nonunique entries')
	return
end

x=x(:); % Turn into column vector

sigmod_result = sigmoid(y(Set).*(x'*X(:,Set)));
sigmod_result = sigmod_result + (sigmod_result<eps).*eps;

F = -log(sigmod_result)/m + lambda*sum(x.^2)/(2*m);

if nargout == 2
    g = -y(Set).*X(:,Set) .* repmat((ones(1,length(Set))-sigmod_result),n,1) ...
        + lambda * repmat(x,1,length(Set));
    J = g/m;
end
end

function [y] = sigmoid(x)
% Sigmoid function.

    y = 1.0 ./ (1.0 + exp(-x));

end