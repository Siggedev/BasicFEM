% VSM167
% Computer Assignment 1
% Written by: Arvid & Linda
clc

% Dimensions

% Widths
w_1 = 1.2; % [m]
w_2 = 0.2; % [m]
w_3 = 0.4; % [m]

% Heights
h_1 = 0.25; % [m]
h_2 = 0.5;  % [m]
h_3 = 0.05; % [m]
h_4 = 0.5;  % [m] 

% Radius
r = 0.1;    % [m] Heat source radius

% Input parameters
T_room = 20; % [C] Room temperature
T_out = -10; % [C] Outside temperature
Q_h = 200; % [W/m] Heat source
q_2 = - Q_h / (2*r*pi); % [W/m2] Heat flux in q_2
nquadpoints = 2; 

% Heat conductivity 
k_floor_tiles = 1.0;     % [W/mC]
k_floor_structure = 2.0; % [W/mC]
k_wall_structure = 0.1;  % [W/mC]

% Heat model
floormesh_tria_path = fullfile("assignments", "CA1", "inputdata", "floormesh_tria.inp");
floormesh_quad_path = fullfile("assignments", "CA1", "inputdata", "floormesh_quad.inp");

% Output from mesh
[element_nodes,     ... % Which nodes belong to each element (connectivity)
node_coordinates,   ... % Where each node is
element_sets,       ... % Groups of elements (regions)
facet_sets,         ... % Groups of edges, as element–facet pairs (boundary pieces)
node_sets,          ... % Groups of nodes (boundary points)
element_type        ... 
] = read_mesh(floormesh_quad_path);


ndofs = size(node_coordinates, 2);

% Ka = f
K = spalloc(ndofs, ndofs, 10*ndofs);
f = zeros(ndofs, 1);
a = zeros(ndofs, 1); 

number_elements = size(element_nodes, 2); 

k_el = zeros(number_elements, 1);
k_el(element_sets{"tiles"}) = k_floor_tiles;
k_el(element_sets{"floorstructure"}) = k_floor_structure;
k_el(element_sets{"wallstructure"}) = k_wall_structure;


for i = 1:number_elements

    element_node_ids = element_nodes(:, i);          % Element node ids
    coords = node_coordinates(:, element_node_ids);  % Element node coordinates

    % Conductivity matrix
    k = k_el(i);
    if k == 0
        error('conductivity cant be zero')
    end
    D = k * eye(2);
    
    % Solving Ke and fe
    switch element_type
        case "LinearTriangle"
            [Ke, fe] = linear_triangle_heat_element(D, 0, coords, nquadpoints);
        case "BilinearQuadrilateral"
            [Ke, fe] = bilinear_quadrilateral_heat_element(D, 0, coords, nquadpoints);
        otherwise
            error("FUCK YOU")
    end

    edofs = element_node_ids;                   % For scalar problems, we use same numbering of nodes and dofs
    K(edofs, edofs) = K(edofs, edofs) + Ke;     % Add stiffness contribution to global matrix
    f(edofs) = f(edofs) + fe;                   % Add load vector contribution to global vector
end

dofs_roomside = union(node_sets{"tilesurface"}, node_sets{"wallsurface"});
dofs_outside = node_sets{"outsurface"};

% Dirichlet boundaries
constrained_dofs = union(dofs_roomside, dofs_outside);

% Neumann boundaries
neumann_dofs = node_sets{"hole"}

% Solve the equation system via partitioning
cdofs = constrained_dofs % Dofs with known temperature, so constrain them!

a(dofs_roomside) = T_room;
a(dofs_outside)  = T_out;

hole_facets = facet_sets{"hole"};

f_neumann = zeros(ndofs, 1);

for i = 1:size(hole_facets, 2)
    
    element_nr = hole_facets(1, i);
    facet_nr = hole_facets(2, i);

    enods = element_nodes(:, element_nr);
    coords = node_coordinates(:, enods);

    switch element_type
        case "LinearTriangle"
            f_neumann(enods) = f_neumann(enods) + linear_triangle_heat_neumann(q_2, coords, facet_nr, nquadpoints);
        case "BilinearQuadrilateral"
            f_neumann(enods) = f_neumann(enods) + bilinear_quadrilateral_heat_neumann(q_2, coords, facet_nr, nquadpoints);
        otherwise
            error("FUCK YOU")
    end


end

% 
% f(neumann_dofs) = q_2;
f = f + f_neumann;
fdofs = setdiff(1:ndofs, cdofs);    % The rest are unknown -> free to change

% Solve for the unknown dofs: af = Kff^{-1} * (ff - Kfc * ac)
af = K(fdofs, fdofs) \ (f(fdofs) - K(fdofs, cdofs) * a(constrained_dofs))

% Re-create the full solution vector
a = zeros(ndofs, 1);
a(fdofs) = af;
a(cdofs) = a(constrained_dofs);