function [residuals] = damped_oscillator_residual(x, Set, t, experimental_values)

    Fx = damped_oscillator(x, Set, t); 

    residuals = Fx - experimental_values(Set); 

end