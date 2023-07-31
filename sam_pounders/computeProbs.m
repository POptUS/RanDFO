function [pi, V] = computeProbs(s,delta)

p = length(delta);

delta = max(delta,eps);

[sorted_delta, sort_inds] = sort(abs(delta)); 
cum_sorted_delta = cumsum(sorted_delta);

% find largest c satisfying lhs <= rhs
lhs = (1:p) + s - p;
rhs = cum_sorted_delta./sorted_delta;
c = find((lhs<=rhs).*(lhs>0),1,'last');

% now compute the probabilities
pi = ones(p,1);
if ~isempty(c)
    pi(sort_inds(1:c)) = (lhs(c)/cum_sorted_delta(c))*sorted_delta(1:c);
end


% to avoid divide by zero issues
pi(pi <= eps) = eps;

V = sum((delta.^2).*(1./pi' - 1.0));

end