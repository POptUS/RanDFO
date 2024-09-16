function z = noisy_rosenbrock(x)
    n = length(x);
    y = 0;
    ef = 0.01;
    for i = 1:(n-1)
        y = y + 100 * (x(i+1) - x(i)^2)^2 + (1 - x(i))^2;
    end
    z = y + ef * (2 * rand - 1);
end