function [to_update, probs, V] = choose_subset(X, Cres, xkin, centers, old_delta, delta, combinemodels, Lip, sense, sketchsize, ng)

    m = length(Cres);
    errors = zeros(1,m);
    x = X(xkin, :);

    gamma = 0.0; % weighting parameter with uniform distribution
    %Lip = max(Lip, delta);

    for i = 1:m
        d = norm(x - X(centers(i), :));
        if strcmp(sense,'m')
            if d == 0 && old_delta(i) == delta
                errors(i) = 0;
            else
                %errors(i) = Lip(i) * (1 + sqrt(n)* Pars(2)) + ... 
                errors(i) = Lip(i) * (d + delta)^2;
                if strcmp(functions(combinemodels).function,'leastsquares')
                    errors(i) = abs(Cres(i)) * errors(i);
                end
                %errors(i) = Lip(i) * (d + delta);
            end
        elseif strcmp(sense,'f')
            if d == 0
                errors(i) = 0;
            else
                %errors(i) = Lip(i) * (1 + sqrt(n)* Pars(2)) ... 
                errors(i) = Lip(i) * d^2;
                if strcmp(functions(combinemodels).function,'leastsquares')
                    errors(i) = abs(Cres(i)) * errors(i);
                end
                %errors(i) = Lip(i) * d;
            end
        end
    end

    if isempty(sketchsize)
%             sortederrors = sort(errors);
%             cumsumerrors = cumsum(sortederrors);
%             mu = 0.5;
%             %c = 0.1;
%             %threshold = mu * c^2 * (delta * ng)^2;
%             threshold = mu * cumsumerrors(m);
%             if ~isinf(threshold)
%                 sketchsize = m - find(cumsumerrors >= threshold, 1, 'first') + 1;
%             else
%                 sketchsize = m;
%             end
%             if isempty(sketchsize)
%                 sketchsize = 1;
%             end
%             sketchsize = min(m, max(sketchsize, 1));
%              [probs, V] = computeProbs(sketchsize, errors);

             mu = 0.05;
             c = 0.0;
             %c = n;
             %c = m * n * min(Lip);
             %kappa = mu * (min(1.0,delta))^4; 
             if strcmp(sense,'f')
                sketchsize = m;
                probs = ones(m,1);
                V = 0;
                %kappa = mu * c^2 * (min(1.0,delta))^4;
                %[sketchsize,probs,V] = chooseSketchSizeandProbs(kappa, errors);
             elseif strcmp(sense,'m')
                kappa = mu * (min(1.0,delta))^4;
                [sketchsize,probs,V] = chooseSketchSizeandProbs(kappa, errors);
             end             
    else       
            %if strcmp(sense, 'm')
            [probs, V] = computeProbs(sketchsize, errors);
            %else
            %sketchsize = m;
            %probs = ones(m,1);
            %V = 0;
            %end        
    end

    % weight probs with uniform distribution?
    uniform = (sketchsize/m)*ones(m,1);
    probs = gamma * uniform + (1.0 - gamma) * probs;

    to_update = realize_random_subset(probs);

    if ~any(isinf(errors))
        new_probs = get_cps_probs(probs', sketchsize);
        if abs(sum(new_probs) - sketchsize) > 1.0
            % THIS IS A SAFEGUARD IF SOMETHING WENT NUMERICALLY AWRY. IT
            % SHOULD TRIGGER RARELY.  
            % 1.0 IS ARBITRARY, BUT INDICATES SOMETHING HAS GONE QUITE
            % WRONG. 
            warning('Something went far too wrong in the CPS computation. Doing a greedy thing on this iteration.')
            [~, to_update] = maxk(new_probs, sketchsize); 
            rejected = false;
        else
            rejected = true;
        end
        while rejected
            % POISSON SAMPLE
            to_update = [];
            for k = 1:m
                % flip coin
                if rand() < new_probs(k)
                    to_update = cat(1,to_update,k);
                end
            end
            if length(to_update) == sketchsize
                rejected = false;
            end
        end
        to_update = to_update(:);
    end
end