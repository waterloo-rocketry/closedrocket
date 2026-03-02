function [matrix_exp] = math_matrix_exp(matrix)
    % matrix exponential series
    matrix_exp = eye(length(matrix)) + matrix + 1/2*matrix^2 + 1/6*matrix^3 + 1/24*matrix^4;% + 1/120*matrix^5;
end