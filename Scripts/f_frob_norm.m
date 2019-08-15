function I = f_frob_norm(matrix_test, matrix_retest)

difference = matrix_test - matrix_retest;
I = vecnorm(difference);
end

