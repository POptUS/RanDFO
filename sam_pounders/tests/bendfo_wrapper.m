function F = bendfo_wrapper(m, n, nprob, x, Set)

    %Set defines the component functions of F you would like 
    if nargin<5 % Use all components
	    Set = 1:length(x); %(in this case there are dim n components)
    elseif length(Set) > length(unique(Set)) 
	    disp('Warning: Set has nonunique entries')
	    return
    end

    if nargin == 5
        F = dfovec(m, n, x, nprob)';

        % trim down - wasteful, but i'm not going to implement Fvec to handle
        % just Set. 
        F = F(Set);
    else
        F = dfovec(m, n, x, nprob);        
    end

end