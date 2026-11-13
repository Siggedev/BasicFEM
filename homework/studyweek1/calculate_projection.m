function [K, f] = calculate_projection(g, shape_function, nquadpoints)
    
    [weights, points] = line_quadrature(nquadpoints);
    
    n = numel(shape_function([0]));

    K = zeros(n,n);
    f = zeros(n,1);

    for q = 1:nquadpoints
        xi_q = points(q);
        w_q = weights(q);
        

        N = shape_function(xi_q);

        K = K + N' * N * w_q;
        f = f + N' * g(xi_q) * w_q;
    end
end