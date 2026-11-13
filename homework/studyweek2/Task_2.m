%% PATHS
clc
here = fileparts(mfilename("fullpath"));

if ~isfolder(fullfile(here, "..", "..", "src"))
    here = fileparts(matlab.desktop.editor.getActiveFilename);
end

addpath(here)
addpath(genpath(fullfile(here, "..", "..", "src")))


%% linear two noded
clc
h = @(x) 1;
x_alpha = [1, 4];
nquadpoints = 2;
shape_function = @linear_line_reference_shape_values;
shape_gradient = @linear_line_reference_shape_gradients;
calculate_element(h, x_alpha, nquadpoints, shape_function, shape_gradient)


%% qudratic 3 noded
clc
h = @(x) 1;
x_alpha = [1, 4, 2.5];
nquadpoints = 3;
shape_function = @quadratic_line_reference_shape_values;
shape_gradient = @quadratic_line_reference_shape_gradients;
calculate_element(h, x_alpha, nquadpoints, shape_function, shape_gradient)

%%  linear two noded
clc

fprintf("Linear two noded:\n\n")

h = @(x) sin(x);
x_alpha = [0, pi];
nquadpoints = [1, 2, 3, 4];
shape_function = @linear_line_reference_shape_values;
shape_gradient = @linear_line_reference_shape_gradients;

for i = 1:numel(nquadpoints)

    I = calculate_element(h, x_alpha, nquadpoints(i), shape_function, shape_gradient);
    fprintf("nquadpoints = %d: \nintegral = %f: \nerror: %f: \n\n", nquadpoints(i), I, I - 2)
end

%%  quadratic three noded
clc

fprintf("quadratic three noded:\n\n")

h = @(x) sin(x);
x_alpha = [0, pi, pi/2];
nquadpoints = [1, 2, 3, 4];
shape_function = @quadratic_line_reference_shape_values;
shape_gradient = @quadratic_line_reference_shape_gradients;

for i = 1:numel(nquadpoints)

    I = calculate_element(h, x_alpha, nquadpoints(i), shape_function, shape_gradient);
    fprintf("nquadpoints = %d: \nintegral = %f: \nerror: %f: \n\n", nquadpoints(i), I, I - 2)

end

%% 2b)
clc
% 1, 2)
k = 1; % thermal conductivity
T_1 = 20;
T_2 = 0;
b = 10e3; % heat source
x_alpha = [0, 0.1, 0.05];
nquadpoints = 2;


% 3)
[Ke, fe] = quadratic_line_heat_element(k, b, x_alpha, nquadpoints);

% 4)

cdofs = [1,2];
fdofs = [3];

a_c = [T_1; T_2];

% ==================================================
K_ff = Ke(fdofs, fdofs);    K_fc = Ke(fdofs, cdofs);
K_cf = Ke(cdofs, fdofs);    K_cc = Ke(cdofs, cdofs);
% ==================================================

f_f = fe(fdofs);

a_f = K_ff \ (f_f - K_fc * a_c);
f_c = K_cf * a_f + K_cc * a_c

% Plotting
figure
a = [a_c; a_f];

xi = linspace(-1, 1, 11);
x = zeros(size(xi));
T = zeros(size(xi));

for i = 1:numel(xi)
    N = quadratic_line_reference_shape_values(xi(i));
    x(i) = N * x_alpha';
    T(i) = N * a;
end

plot(x, T);
xlabel("\xi")
ylabel("Temperature")
title("Distributed temperature")
grid on


%% 2c)

clc
% 1, 2)
k = 1; % thermal conductivity
T_1 = 20;
b = 10e3; % heat source
q_n = 200; % heat flux
x_alpha = [0, 0.1, 0.05];
nquadpoints = 2;

[Ke, fe] = quadratic_line_heat_element(k, b, x_alpha, nquadpoints);

cdofs = [1];
fdofs = [2,3];

a_c = [T_1];

K_ff = Ke(fdofs, fdofs);    K_fc = Ke(fdofs, cdofs);
K_cf = Ke(cdofs, fdofs);    K_cc = Ke(cdofs, cdofs);

fe(2) = fe(2) - q_n;

f_f = fe(fdofs);

a_f = K_ff \ (f_f - K_fc * a_c);
f_c = K_cf * a_f + K_cc * a_c

% Plot
figure
a = [a_c; a_f];

xi = linspace(-1, 1, 11);
x = zeros(size(xi));
T = zeros(size(xi));

for i = 1:numel(xi)
    N = quadratic_line_reference_shape_values(xi(i));
    x(i) = N * x_alpha';
    T(i) = N * a;
end

plot(x, T);
xlabel("x")
ylabel("Temperature")
title("Distributed temperature")
grid on
