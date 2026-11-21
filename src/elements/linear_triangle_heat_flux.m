function q = linear_triangle_heat_flux(D, coords, ae, nquadpoints)
% q = linear_triangle_heat_flux(D, coords, ae, nquadpoints)
% Calculate the flux in each quadrature point for a linear triangle
% Inputs: D [2,2]       - conductivity tensor
%         coords [2,3]  - element coordinates
%         ae [3,1]      - element temperature dofs
%         nquadpoints   - number of quadrature points
% Output: q [2, nquadpoints] - flux in each quadrature point
    
    [~, points]= triangle_quadrature(nquadpoints);

    q = zeros(2, nquadpoints);   

    for i = 1:nquadpoints
        
        xi = points(:,i);
        
        dN_dxi = linear_triangle_reference_shape_gradients(xi);
        
        J = calculate_jacobian(dN_dxi, coords);
        
        q(:, i ) = - D * (J' \ dN_dxi) * ae;

    end
end
