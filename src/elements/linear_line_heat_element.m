function [Ke, fe] = linear_line_heat_element(D, b, coords, nquadpoints)
% [Ke, fe] = linear_line_heat_element(D, b, coords, nquadpoints)
% 
% Calculate the element stiffness and load vector for heat conduction 
% on a linear line element, with constant conductivity, `D`, and 
% constant heat (per length), `b`. Use `nquadpoints` to perform the
% numerical quadrature.
% Inputs
% * D: [1,1] - conductivity following Fourier's law
% * b: [1] - heat source per length
% * coords: [1,2] - nodal coordinates of the element
% * nquadpoints: [1] - Number of quadrature points for numerical
%   integration
% Outputs
% * Ke: [2,2] - the element stiffness matrix
% * fe: [2,1] - the element load vector (due to heat source)
    
    Ke = zeros(2, 2);
    fe = zeros(2, 1);

    [w_q, xi_q] = line_quadrature(nquadpoints);

    for i = 1:nquadpoints
        
        N = linear_line_reference_shape_values(xi_q(i));
    
        dNdxi = linear_line_reference_shape_gradients(xi_q(i));

        J = calculate_jacobian(dNdxi, coords);

        Ke = Ke + dNdxi' * inv(J) * D * dNdxi * inv(J) * J * w_q(i);
        fe = fe + N' * b * J * w_q(i);

    end

end
