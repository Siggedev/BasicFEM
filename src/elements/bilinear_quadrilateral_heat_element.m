function [Ke, fe] = bilinear_quadrilateral_heat_element(D, b, coords, nquadpoints)
% [Ke, fe] = bilinear_quadrilateral_heat_element(D, b, coords, nquadpoints)
% 
% Calculate the element stiffness and load vector for heat conduction 
% on a bilinear quadrilateral element, with constant conductivity, `D`, and 
% constant heat (per area), `b`. Use `nquadpoints` to perform the
% numerical quadrature.
% Inputs
% * D: [2,2] - conductivity following Fourier's law
% * b: [1] - heat source per length
% * coords: [2,4] - nodal coordinates of the element
% * nquadpoints: [1] - Number of quadrature points in each direction for 
%   numerical integration
% Outputs
% * Ke: [4,4] - the element stiffness matrix
% * fe: [4,1] - the element load vector (due to heat source)
    
    [w_q, xi_q] = quadrilateral_quadrature(nquadpoints);


    Ke = zeros(4, 4);
    fe = zeros(4, 1);

    for i = 1:nquadpoints
        
        xi = xi_q(i);
        w = w_q(i);

        N = bilinear_quadrilateral_reference_shape_values(xi);

        dNdxi = bilinear_quadrilateral_reference_shape_gradients(xi);
        
        J = 


    
    
    end




end
