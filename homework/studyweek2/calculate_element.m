function L = calculate_element( ...
    h, ...
    x_alpha, ...
    nqp, ...
    shape_function, ...
    shape_gradients)
    
    [w_q, xi_q] = line_quadrature(nqp); % Get points and weights
    
    L = 0; % Length of element

    for i = 1:numel(xi_q)

        N = shape_function(xi_q(i)); % Shape value

        dNdxi = shape_gradients(xi_q(i)); % Gradient value
        
        x_q = N * x_alpha'; % Physical coordinate x along the element
        
        J = calculate_jacobian(dNdxi, x_alpha);

        L = L + h(x_q) * J * w_q(i);
        
    end
    
end