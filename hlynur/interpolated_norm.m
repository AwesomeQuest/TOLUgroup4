function n = interpolated_norm(x, y)
    % Interpolated Norm Function
    % This function takes two vectors x and y of different lengths,
    % interpolates them onto a common grid, and computes the norm of the difference.

    % Define the sample positions for x and y
    xi = linspace(0, 1, length(x));
    yi = linspace(0, 1, length(y));

    % Define the common grid t
    N = max(length(x), length(y));
    t = linspace(0, 1, N);

    % Interpolate x and y onto t
    x_interp = interp1(xi, x, t, 'linear');
    y_interp = interp1(yi, y, t, 'linear');

    % Compute the norm of the difference
    n = norm(x_interp - y_interp);
end
