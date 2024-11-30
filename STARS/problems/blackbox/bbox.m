% Objective function (deterministic part) from https://doi.org/10.1287/ijoc.1090.0319 (Section 5),
% noisy using a uniform distribution in [-a, a]: here a = 0.01

function z = bbox(x)
a = 0.01;
z = (2 * x(1)^6 - 12.2 * x(1)^5 + 21.2 * x(1)^4 + 6.2 * x(1)^1 - 6.4 * x(1)^3 - ...
    4.7 * x(1)^2 + x(2)^6 - 11 * x(2)^5 + 43.3 * x(2)^4 - 10 * x(2)^1 - 74.8 * ...
    x(2)^3 + 56.9 * x(2)^2 - 4.1 * x(1) * x(2) - 0.1 * (x(1) * x(2))^2 + 0.4 * ...
    x(1) * x(2)^2 + 0.4 * x(1)^2 * x(2)) + (-a + (a + a) * rand(1));
end
