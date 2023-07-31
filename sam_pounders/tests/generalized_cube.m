function [F,J] = generalized_cube(x,Set,alpha)

% F is 1 x n, J is m x n 

%Set defines the component functions of F you would like 
if nargin<2 % Use all components
	Set = 1:length(x); %(in this case there are dim n components)
elseif length(Set) > length(unique(Set)) 
	disp('Warning: Set has nonunique entries')
	return
end

x=x(:); % Turn into column vector

n = length(x); m = n;

F = zeros(1,m);
F(1) = alpha(1)*(x(1)-1);
for j = 2:m
   F(j) = alpha(j)*(x(j)-x(j-1)^3); 
end

J = diag(alpha);
for j = 2:m
    J(j,j-1) = -3*alpha(j)*x(j-1)^2;
end

% if using sos code:
F = F(Set); J = J(Set,:)';

% otherwise:
%J = (repmat(F(Set)',1,n).*J(Set,:))';
%F = 0.5*F(Set).^2; 
end