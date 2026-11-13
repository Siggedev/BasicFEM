%% Paths
% ==============================================
% RUN THIS SECTION BEFORE RUNNING OTHER SECTIONS
% ==============================================


% Folder of this script
here = fileparts(mfilename('fullpath'));

% Run Section executes a temporary copy -> fall back to the file open in the Editor
if ~isfolder(fullfile(here, '..', '..', 'src'))
    here = fileparts(matlab.desktop.editor.getActiveFilename);
end

% Add this script's folder (local functions) and BasicFEM/src (with subfolders) to the path
addpath(here)
addpath(genpath(fullfile(here, '..', '..', 'src')))

%%
%  1a)
clc
figure
hold on

xi = linspace(-1, 1, 101); % 101 => Even decimal number 

a_1 = [1, 0, 0]';
a_2 = [0, 1, 0]';
a_3 = [0, 0, 1]';
a_4 = [0.5, 1, 0.75]';
a_5 = [1, 1, 1]';
a_6 = [0, 1, 1]';

a = [a_1, a_2, a_3, a_4, a_5, a_6];

ncols = size(a, 2);

for j = 1:ncols
    
    f = zeros(size(xi));
    aj = a(:,j);

    for i = 1:numel(xi)
        N    = quadratic_line_reference_shape_values(xi(i));
        f(i) = N(1) * aj(1) + N(2) * aj(2) + N(3) * aj(3);
    end

    plot(xi, f)
    grid on
    xlabel('\xi')
    ylabel('f(\xi)')

end
hold off

%% 1b)
clc

% Domain
xi = linspace(-1, 1, 101);

% Define figure
figure

% Hold plotting
hold on

% Define function and plot true values
g = @(xi) exp(xi);
plot(xi, g(xi), "--", LineWidth=2);

% Iterate over quadrature points 3 and 4
for i = 3:4
    
    % calculate K and f matrix
    [K, f] = numerical_integration(i);

    % Find coefficients 
    a = K \ f;
    
    % Approximated values
    % Placeholder array with zeros
    u = zeros(size(xi));

    % Iterate over domain values
    for j = 1:numel(xi)

        % Calculate shape values
        N = quadratic_line_reference_shape_values(xi(j));

        % Add approximated value 
        u(j) = N(1)*a(1) + N(2)*a(2) + N(3)*a(3);
    end
    
    % Plot approximated function
    plot(xi, u);

end

% Plot functions
xlabel("\xi")
grid on 
legend("e^\xi", "3 quad points", "4 quad points")
hold off

%% 1c)
clc

% DOMAIN
xi = linspace(-1, 1, 101);

% INPUTS
g = @(xi) exp(xi);
shape_functions = {@quadratic_line_reference_shape_values, @linear_line_reference_shape_values};
nquadpoints = 4;

figure
hold on
plot(xi, g(xi), "--", LineWidth=2)

% K and f matrix

for k = 1:numel(shape_functions)
    shape_func = shape_functions{k};
    [K, f] = calculate_projection(g, shape_func, nquadpoints);
    a = K \ f;

    u = zeros(size(xi));

    for i = 1:numel(xi)

        N = shape_func(xi(i));

        u(i) = N * a;
    end

    plot(xi, u)
end

xlabel("\xi")
grid on
legend("e^\xi", "quadratic line function", "linear line function")
hold off
