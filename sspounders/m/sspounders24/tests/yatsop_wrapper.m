function F = yatsop_wrapper(x,probspecs,probtype)

    [~,F] = calfun_sample(x', probspecs, probtype);
    F = F';

end