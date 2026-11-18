
path = genpath("src");
addpath(path)

%%
clc

x1 = 1; x2 = 1.5; x3 = 2;
e = 1;

node_coordinates = [x1, x2, x3];

e1_nodes = [1; 2]
e2_nodes = [2; 3]
element_nodes = [e1_nodes, e2_nodes]

enods = element_nodes(:, e)
ecoords = node_coordinates(:, enods)


%% 

nel = size(element_nodes, 2);        % Number of elements
ndofs = size(node_coordinates, 2);   % Number of nodes
K = spalloc(ndofs, ndofs, 3 * ndofs);% Allocate global matrix
f = zeros(ndofs, 1);				 % Allocate global vector

for e = 1:nel
    enods = element_nodes(:, e);     % Extract node numbers for given element
    ecoords = node_coordinates(:, enods); % Extract node coordinates for elem
    % Calculate the local matrix and vector for the given element:
    [Ke, fe] = linear_line_heat_element(k, b, ecoords, nquadpoints);
    % Add the local contributions to the right locations in the global 
    % matrix and vector:
    edofs = enods;
    K(edofs, edofs) = K(edofs, edofs) + Ke;
    f(edofs) = f(edofs) + fe;
end