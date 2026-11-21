function fe = bilinear_quadrilateral_heat_neumann(qn, coords, facet_nr, nquadpoints)
% function fe = bilinear_quadrilateral_heat_neumann(qn, coords, facetnr, nquadpoints)
% Calculate the contribution to the element load vector due to a constant 
% outwards boundary flux, `qn`, on the facet, `facetnr`, of a bilinear 
% quadrilateral element.
% Input: qn          - outwards heat flux
%        coords[2x4] - element coordinates
%        facet_nr     - The facet number where the heat flux is applied
%        nquadpoints - Number of quadrature points to use for numerical
%                      integration
% Output: fe [4x1]   - Local element load vector contribution
    [weights, points] = line_quadrature(nquadpoints);
    
    fe = zeros(4, 1);
    
    for i = 1:nquadpoints
    
        xi_1D = points(:,i);
        w = weights(i);
    
        xi_2D = quadrilateral_facet_coords(xi_1D, facet_nr);
    
        N = bilinear_quadrilateral_reference_shape_values(xi_2D);
        dN_dxi = bilinear_quadrilateral_reference_shape_gradients(xi_2D);
    
        J = calculate_jacobian(dN_dxi, coords);
    
        n = quadrilateral_facet_weighted_normal(J, facet_nr);
    
        fe = fe - N' * qn * norm(n) * w;
    end
end


