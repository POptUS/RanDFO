function F = calfun_wrapper(x)
    
    global BenDFO

    [~, F] = calfun(x, BenDFO, 'smooth');

    F = F';

end