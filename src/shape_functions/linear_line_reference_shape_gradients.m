function dNdxi = linear_line_reference_shape_gradients(xi)
% dNdxi = linear_line_reference_shape_gradients(xi)
% Calculate the shape function gradients for a linear line element. 
% Input: xi [1x1] - local coordinate
% Output: dNdxi [1x2] - shape function gradients at xi
%         dNdxi(:, i) is the gradient of the ith shape function
    if ~all(size(xi) == 1)
        error('xi must be 1x1 matrix')
    end

    dNdxi = [-1/2, 1/2];
end

