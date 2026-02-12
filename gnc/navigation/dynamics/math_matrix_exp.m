function [matrix_exp] = math_matrix_exp(matrix)
    % matrix exponential series without unity term
    matrix_exp = matrix + 1/2*matrix^2 + 1/6*matrix^3 + 1/24*matrix^4;
end