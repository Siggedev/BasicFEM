function equal = my_isapprox(a, b, rtol, atol)
    if nargin < 4 || isempty(atol)
        atol = 0;
    end
    if nargin < 3 || isempty(rtol)
        if atol > 0
            rtol = 0;
        else
            rtol = sqrt(eps);
        end
    end

    if ~isequal(size(a), size(b))
        equal = false;
        return
    end

    % Flatten to ensure Frobenius norm for matrices
    na = norm(a(:));
    nb = norm(b(:));
    nd = norm(a(:) - b(:));

    if isfinite(nd)
        equal = nd <= max(atol, rtol * max(na, nb));
    else
        % Fallback to element-wise comparison
        tol = max(atol, rtol .* max(abs(a), abs(b)));
        equal = (a == b) | (abs(a - b) <= tol);
        equal = all(equal(:));
    end
end
