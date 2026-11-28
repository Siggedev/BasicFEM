function q = bilinear_quadrilateral_heat_flux(D, coords, ae, nquadpoints)
% q = bilinear_quadrilateral_heat_flux(D, coords, ae, nquadpoints)
% Calculate the flux in each quadrature point for a bilinear quadrilateral
% Inputs: D [2,2]       - conductivity tensor
%         coords [2,4]  - element coordinates
%         ae [4,1]      - element temperature dofs
%         nquadpoints   - number of quadrature points in each direction
% Output: q [2, nquadpoints^2] - flux in each quadrature point
    [~, points]= quadrilateral_quadrature(nquadpoints);
    
    q = zeros(2, nquadpoints^2);   
    
    for i = 1:nquadpoints^2
    
        xi = points(:,i);
    
        dN_dxi = bilinear_quadrilateral_reference_shape_gradients(xi);
    
        J = calculate_jacobian(dN_dxi, coords);
    
        q(:, i) = - D * (J' \ dN_dxi) * ae;
    
    end
end
