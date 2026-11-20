function [Ke, fe] = linear_triangle_heat_element(D, b, coords, nquadpoints)
% [Ke, fe] = linear_triangle_heat_element(D, b, coords, nquadpoints)
% 
% Calculate the element stiffness and load vector for heat conduction 
% on a linear triangle element, with constant conductivity, `D`, and 
% constant heat (per area), `b`. Use `nquadpoints` to perform the
% numerical quadrature.
% Inputs
% * D: [2,2] - conductivity following Fourier's law
% * b: [1] - heat source per length
% * coords: [2,3] - nodal coordinates of the element
% * nquadpoints: [1] - Number of quadrature points for numerical
%   integration
% Outputs
% * Ke: [3,3] - the element stiffness matrix
% * fe: [3,1] - the element load vector (due to heat source)
    
    [weights, points] = triangle_quadrature(nquadpoints);

    Ke = zeros(3,3);
    fe = zeros(3,1);

    for i = 1:nquadpoints
        
        xi = points(:,i);
        w = weights(i);

        N = linear_triangle_reference_shape_values(xi);

        dNdxi = linear_triangle_reference_shape_gradients(xi);

        J = calculate_jacobian(dNdxi, coords);
            
        B = J' \ dNdxi;
        
        Ke = Ke + B' * D * B * det(J) * w;

        fe = fe + N' * b * det(J) * w;
    end
end
