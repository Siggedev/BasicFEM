function [weights, points] = quadrilateral_quadrature(nquadpoints)
% [weights, points] = quadrilateral_quadrature(nquadpoints)
% Return the weights and quadrature points for numerically integrating
% the reference quadrilateral with `nquadpoints^2` points.
% Input: nquadpoints [1x1] Number of quadrature points in each direction
% Output: weights [1, nquadpoints^2] - quadrature weights
%         points  [2, nquadpoints^2] - quadrature points
    
    % Quadrilateral rule:
    % 1D Gauss rule applied in both directions (a tensor product)
    
    [w_1d, p_1d] = line_quadrature(nquadpoints);

    n = nquadpoints^2;
    
    weights = zeros(1, n);
    points  = zeros(2, n);
    
    q = 1;

    for i = 1:nquadpoints

        for j = 1:nquadpoints
            
            weights(q)   = w_1d(j) * w_1d(i);
            points(:,q)  = [p_1d(j); p_1d(i)];

            q = q + 1;
        end     
    
    end

end
    

