function [mval, grad, Hess] = MinFrobQuadMod(Y, F)
    % Y - a m-by-n matrix of m data points of dimension n: Y = [y1; y2; ...; ym]
    % that is, yi is an n-dimensional vector, for i = 1, 2, ..., m
    % F is an m-dimensional vector containing the values of the function f
    % at the m interpolation points yi, i = 1, 2, ..., m
    
    if isrow(F)
        F = F';
    end
    
    [m, n] = size(Y);
    
    K = zeros((n+1)*(n+2)/2, (n+1)*(n+2)/2);
    
    K((n+2):end, (n+2):end) = eye(n*(n+1)/2);
    
    B = [ones(m, 1), Y, (1/2) * matYYT(Y)];
    
    zvec = zeros((n+1)*(n+2)/2, 1);

    options = optimoptions('quadprog', 'Display', 'off');

    x = quadprog(K, zvec, [], [], B, F, [], [], [], options);
    
    mval = x(1);
    
    grad = x(2:(n+1));
    
    h = x((n+2):end);
    
    Hess = vecToSymmetric(h, n);
end

