function [K, f] = numerical_integration(nquadpoints)
   % function [K, f] = numerical_integration(nquadpoints)
   % Performs numerical integration
   % Input: nquadpoints [1x1] Number of quadrature points
   % Output: calculated K & f matrix
   
   [weights, points] = line_quadrature(nquadpoints);
    
   K = spalloc(3, 3, 9);
   f = spalloc(3, 1, 3); 

   for q = 1:nquadpoints
        xi_q = points(q);
        w_q = weights(q);

        N = quadratic_line_reference_shape_values(xi_q);

        K = K + N' * N * w_q;
        f = f + N' * exp(xi_q) * w_q;
   end
end