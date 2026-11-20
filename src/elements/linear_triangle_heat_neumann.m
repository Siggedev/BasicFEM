function fe = linear_triangle_heat_neumann(qn, coords, facet_nr, nquadpoints)
% function fe = linear_triangle_heat_neumann(qn, coords, facetnr, nquadpoints)
% Calculate the contribution to the element load vector due to a constant 
% outwards boundary flux, `qn`, on the facet, `facet_nr`, of a linear 
% triangle element.
% Input: qn          - outwards heat flux
%        coords[2x3] - element coordinates
%        facet_nr     - The facet number where the heat flux is applied
%        nquadpoints - Number of quadrature points to use for numerical
%                      integration
% Output: fe [3x1]   - Local element load vector contribution
    
% triangle_facet_coords(xi_line, face_nr)
% triangle_facet_weighted_normal(J, facet_nr)

    [weights, points] = line_quadrature(nquadpoints);

    fe = zeros(3, 1);
    
    for i = 1:nquadpoints
        
        xi_1D = points(:,i);
        w = weights(i);

        xi_2D = triangle_facet_coords(xi_1D, facet_nr);
        
        N = linear_triangle_reference_shape_values(xi_2D);
        dN_dxi = linear_triangle_reference_shape_gradients(xi_2D);
        
        J = calculate_jacobian(dN_dxi, coords);
        
        n = triangle_facet_weighted_normal(J, facet_nr);

        fe = fe - N' * qn * norm(n) * w;
    end
end
