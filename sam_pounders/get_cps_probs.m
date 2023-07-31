function [pi_tilde,lambda,jip] = get_cps_probs(pi,b)

% COMPUTES DISTRIBUTION FOR CONDITIONAL POISSON SAMPLING THAT MATCHES THE
% INCLUSION PROBABILITIES FOR POISSON SAMPLING 

% fixing to avoid numerical frustrations:
pi(pi >= 1.0 - sqrt(eps)) = 1.0 - sqrt(eps);

% initial point
pi_tilde = pi; p = length(pi_tilde);

residual_old = inf; maxiters = 1000; iter = 0;
psi = recursive(pi_tilde,b);
while residual_old > p*sqrt(eps) && iter < maxiters    
    alpha = 1.0; residual_trial = inf;
    %jip = grad_recursive(pi_tilde,psi);
    %Delta = jip - (pi_tilde'*pi_tilde);
    %J = Delta*diag(1.0./(pi_tilde.*(1.0-pi_tilde)));
    % USING AN IDENTITY APPROXIMATION FOR J:
    while residual_trial >= residual_old && alpha > eps
        pi_tilde_new = pi_tilde + alpha*(pi-psi);
        %pi_tilde_new = pi_tilde + alpha*(J\(pi-psi)')';
        if all(pi_tilde_new > 0) && all(pi_tilde_new < 1)
            psi_trial = recursive(pi_tilde_new,b);
            residual_trial = norm(pi-psi_trial);
        end
        alpha = 0.5*alpha;
    end
    pi_tilde = pi_tilde_new;

    psi = psi_trial;
    residual_old = norm(pi-psi);
    iter = iter + 1;
end

%% OUTPUT CHECKING
% is pi_tilde a proper distribution? 
if any(pi_tilde <= 0) || any(pi_tilde >= 1) || any(isnan(pi_tilde)) || any(~isreal(pi_tilde))
    % unfortunate! 
    warning('There was an error in updating pi to pi_tilde.')
    pi_tilde = pi;
end

% now pick a scalar c so that pi_tilde sums to b: 
lambda = log(pi_tilde./(1.0-pi_tilde));

% solve sum(exp(lambda + c)./(1.0 + exp(lambda+c))) = b by Newton method;
fun = @(c)sumexplambda(c,lambda);
residual_old = inf; maxiters = 1000; iter = 0; c = 0;
[fc,gc] = fun(c); 
while residual_old > sqrt(eps) && iter < maxiters    
    alpha = 1.0;  residual_trial = inf;
    while residual_trial >= residual_old && alpha > eps
        cnew = c - alpha*(fc-b)/gc;
        [fcnew,gcnew] = fun(cnew);
        residual_trial = abs(fcnew-b);
        alpha = 0.5*alpha;
    end
    c = cnew;
    fc = fcnew; gc = gcnew;
    residual_old = abs(fc-b);
    iter = iter + 1;
end

[~,~,pi_tilde] = fun(c);

%% OUTPUT CHECKING
% is pi_tilde a proper distribution? 
if any(pi_tilde < 0) || any(pi_tilde > 1) || any(isnan(pi_tilde)) || any(~isreal(pi_tilde))
    warning('There was an error in updating pi to pi_tilde.')
    pi_tilde = pi;
end
end % end main function

function [psi] = recursive(pi_tilde,b)
    lower_psi = zeros(1,length(pi_tilde));
    scaling = pi_tilde./(1.0-pi_tilde);
    for k = 1:b    
        psi_num = scaling.*(1.0 - lower_psi);
        psi_denom = sum(psi_num);
        lower_psi = (1/psi_denom)*psi_num;
        if any(lower_psi >= 1)           
            lower_psi(lower_psi >= 1) = 1 - eps;
        end
    end
    psi = b*lower_psi;
end

function [jip] = grad_recursive(pi_tilde,psi)
    lambda = log(pi_tilde./(1.0-pi_tilde));
    p = length(lambda);
    jip = zeros(p); dolater_k = []; dolater_l = []; 
    for k = 1:p
        for l = (k+1):p
            if abs(lambda(k) - lambda(l)) > sqrt(eps)
                jip(k,l) = (psi(k)*exp(lambda(l)) - psi(l)*exp(lambda(k)))/(exp(lambda(l))-exp(lambda(k)));
            else
                dolater_k = cat(1,dolater_k,k); dolater_l = cat(1,dolater_l,l);
            end
        end
    end
    jip = jip + jip'; jip = jip + diag(psi);
    remaining = length(dolater_k);
    if remaining > 0
        for r = 1:remaining
            k = dolater_k(r); l = dolater_l(r);
            inds = find(abs(lambda - lambda(l)) <= sqrt(eps));
            skippedinds = setdiff(1:p,inds);
            sumnoneq = sum(jip(k,skippedinds));
            jip(k,l) = 1/(length(inds)-1)*(psi(l)*(b-1) - sumnoneq);
            jip(l,k) = jip(k,l);
        end
    end
end

function [f,g,pi_tilde] = sumexplambda(c,lambda)
    u = exp(lambda + c);
    f = sum(u./(1.0 + u));
    g = sum((u.*(1.0 + u) - u.^2)./((1.0 + u).^2));
    pi_tilde = u./(1.0 + u);
end