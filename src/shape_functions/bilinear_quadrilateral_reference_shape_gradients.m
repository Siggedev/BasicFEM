function dNdxi = bilinear_quadrilateral_reference_shape_gradients(xi)
% dNdxi = bilinear_quadrilateral_reference_shape_gradients(xi)
% Calculate the shape function gradients for a bilinear quadrilateral elem.
% Input: xi [2x1] - local coordinates
% Output: dNdxi [2x4] - shape function gradients at xi
%         dNdxi(:, i) is the gradient of the ith shape function
    if size(xi,1) ~= 2 || size(xi,2) ~= 1
        error("xi must be 2x1")
    end
    
    xi_1 = xi(1);
    xi_2 = xi(2);

    dNdxi_1_1 = -1 + xi_2;
    dNdxi_1_2 = 1 - xi_2;
    dNdxi_1_3 = 1 + xi_2;
    dNdxi_1_4 = -1 - xi_2;

    dNdxi_2_1 = -1 + xi_1;
    dNdxi_2_2 = -1 - xi_1;
    dNdxi_2_3 = 1 + xi_1;
    dNdxi_2_4 = 1 - xi_1;

    dNdxi = 0.25 * [dNdxi_1_1, dNdxi_1_2, dNdxi_1_3, dNdxi_1_4;
                    dNdxi_2_1, dNdxi_2_2, dNdxi_2_3, dNdxi_2_4];
end
