function N = linear_triangle_reference_shape_values(xi)
% N = linear_triangle_reference_shape_values(xi)
% Calculate the shape function values for a linear triangle element. 
% Input: xi [2x1] - local coordinates
% Output: N [1x3] - shape function values at xi
    if size(xi,1) ~= 2 || size(xi,2) ~= 1
        error("xi must be 2x1")
    end
    N = [xi(1), xi(2), 1 - xi(1) - xi(2)];
end
