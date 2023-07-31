function [sk,pi,V] = chooseSketchSizeandProbs(kappa,errors)
    
    p = length(errors);

    if p == 0
        sk = 0;
        pi = [];
        V = Inf;
    else

        b = 1;
        while b <= p 
            pi = computeProbs(b,errors);
            V = sum((errors.^2).*(1./pi' - 1.0));
            if V <= kappa
                break;
            else   
                b = b + 1;
            end
        end
    
        sk = b;

    end

end