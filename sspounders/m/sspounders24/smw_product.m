function Ainvb = smw_product(sketch_dictionary, regularizer, b)

% implements the SMW rank one updates using the sketches in sketch
% dictionary

% b is what the inverse covariance is being multiplied by

num_sketches = length(sketch_dictionary);

for j = 1:num_sketches
    sketch_dictionary{j}.AinvU = (1.0 / regularizer) * sketch_dictionary{j}.S;
end
    
Ainvb = (1.0 / regularizer) * b;

for i = 1:num_sketches
    U = sketch_dictionary{i}.S;
    U_size = sketch_dictionary{i}.sketch_size;
    AinvU = sketch_dictionary{i}.AinvU;

    inner_product = inv(eye(U_size) + U' * AinvU) * U';
    Ainvb = Ainvb - AinvU * inner_product * Ainvb;

    for j = (i + 1):num_sketches
        AinvU_plus = sketch_dictionary{j}.AinvU;
        sketch_dictionary{j}.AinvU = AinvU_plus - AinvU * inner_product * AinvU_plus;
    end
end

end

