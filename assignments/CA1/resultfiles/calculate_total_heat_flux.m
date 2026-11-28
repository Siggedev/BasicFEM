function Q = calculate_total_heat_flux( ...
    facet_sets, ...
    element_nodes, ...
    element_type, ...
    node_coordinates, ...
    nquadpoints, ...
    q, ...
    boundary_segment)
    
    
    % [element_nr, facet_nr] 
    elements = facet_sets{boundary_segment};
    element_numbers = elements(1,:);

    % Flux in [U;V] direction
    q_n = q(:, elements(1,:));
    
    % Get element node numbers
    element_node_numbers = element_nodes(:, element_numbers);
    
    % Number of elements in tilesurface
    num_elements = size(elements, 2);

    switch element_type
        case "LinearTriangle"
            [~, points] = triangle_quadrature(nquadpoints);
        case "BilinearQuadrilateral"
            [~, points] = quadrilateral_quadrature(nquadpoints);
    end

    Q = 0;

    for i = 1:num_elements

        switch element_type
            case "LinearTriangle"
                dNdxi = linear_triangle_reference_shape_gradients(points(:,1));
                coords = node_coordinates(:, element_node_numbers(:,i));
                J = calculate_jacobian(dNdxi, coords);
                facet_nr = elements(2, i);
    
                n_w = triangle_facet_weighted_normal(J, facet_nr);
                Q = Q + 2 * q_n(:,i)' * n_w;
            case "BilinearQuadrilateral"
        end

    end

end