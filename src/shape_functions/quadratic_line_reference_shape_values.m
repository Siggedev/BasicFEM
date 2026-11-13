function N = quadratic_line_reference_shape_values(xi)
% N = linear_triangle_reference_shape_values(xi)
% Calculate the shape function values for a quadratic line element. 
% Input: xi [1x1] - local coordinate
% Output: N [1x3] - shape function values at xi

    if ~all(size(xi) == 1)
        error("xi must be 1x1")
    end

    N_1 = xi(1) * (xi(1) - 1) / 2;
    N_2 = xi(1) * (xi(1) + 1) / 2;
    N_3 = 1 - xi(1)^2;

    N = [N_1, N_2, N_3];

end
