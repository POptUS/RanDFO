function [F,J] = basic_cubic(x,Set,alpha)

% F is 1 x n, J is m x n 

%Set defines the component functions of F you would like 
if nargin<2 % Use all components
	Set = 1:length(x); %(in this case there are dim n components)
elseif length(Set) > length(unique(Set)) 
	disp('Warning: Set has nonunique entries')
	return
end

n = length(x); m = n;

F = x + (1/2.0) * x.^2 + (1/3.0) * (alpha.*(abs(x).^3));

J = ones(m, n) + repmat(x, m, 1) + alpha(:) * (x.^2);

% if using sos code:
F = F(Set); J = J(Set,:)';

end