function [Lip, old_centers, old_gs] = update_lipschitz(Lip, X, xkin, old_centers, centers, old_gs, Gres)

    m = length(Lip);

    for i = 1:m
        oci = old_centers(i);
        if oci == 0 && centers(i) == xkin %valid(i)
            old_centers(i) = xkin;
            old_gs(:, i) = Gres(:, i);
        elseif oci ~= xkin && centers(i) == xkin %valid(i)
            num = norm(Gres(:, i) - old_gs(:, i));
            denom = norm(X(xkin, :) - X(oci, :));
            rho = num/denom;
            %rho = 1.0;
            if ~isinf(rho) && ~isnan(rho)
                if isinf(Lip(i))
                    Lip(i) = abs(rho);
                    %Lip(i) = norm(Gres(:, i));
                else
                    Lip(i) = max(Lip(i), abs(rho)); 
                    %Lip(i) = abs(rho);
                    %Lip(i) = norm(Gres(:, i));
                end
            end
            old_gs(:, i) = Gres(:, i);
            old_centers(i) = xkin;
        else
            % make more likely to be sampled in the future
            %Lip(i) = 1.1 * Lip(i);
        end
    end

end