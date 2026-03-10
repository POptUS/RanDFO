function [F,J] = generalized_rosenbrock(x,Set,alpha)

% F is 1 x m, J is m x n 

%Set defines the component functions of F you would like 
if nargin<2 % Use all components
	Set = 1:length(x); %(in this case there are dim n components)
elseif length(Set) > length(unique(Set)) 
	disp('Warning: Set has nonunique entries')
	return
end

x=x(:); % Turn into column vector

n = length(x); m = 2 * (n - 1);

F = zeros(1, m); J = zeros(m);
for j = 1:(m/2)
   F(j) = 10*alpha(j)*(x(j+1)^2 - x(j));
   J(j,j+1) = 20*alpha(j)*x(j+1);
   J(j,j) = -10*alpha(j);
end
for j = (m/2 + 1):m
    F(j) = alpha(j) * (1-x(j - (m/2)));
    J(j,j) = -alpha(j);
end

% if using sos code:
F = F(Set); J = J(Set,:)';

% otherwise:
%J = (repmat(F(Set)',1,n).*J(Set,:))';
%F = 0.5*F(Set).^2; 
end