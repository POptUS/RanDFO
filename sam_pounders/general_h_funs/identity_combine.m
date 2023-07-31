function [G, H] = identity_combine(Cres, Gres, Hres)

    G = sum(Gres,2);
    H = sum(Hres,3);
    
end