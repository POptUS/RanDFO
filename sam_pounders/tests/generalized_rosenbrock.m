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

n = length(x); m = n;

if mod(n,2) ~= 0
    error('n must be even');
end

F = zeros(1,m); J = zeros(m);
for j = 1:(m/2)
   F(j) = 10*alpha(j)*(x(2*j-1)^2 - x(2*j));
   F((m/2)+j) = alpha((m/2)+j)*(x(2*j-1)-1);
   J(j,2*j-1) = 20*alpha(j)*x(2*j-1);
   J(j,2*j) = -10*alpha(j);
   J((m/2)+j,2*j-1) = alpha((m/2)+j);
end


% if using sos code:
F = F(Set); J = J(Set,:)';

% otherwise:
%J = (repmat(F(Set)',1,n).*J(Set,:))';
%F = 0.5*F(Set).^2; 
end