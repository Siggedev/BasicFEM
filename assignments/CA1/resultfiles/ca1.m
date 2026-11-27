% VSM167
% Computer Assignment 1
% Written by: Arvid & Linda

% Add src
addpath(genpath("src"))

%% Solve FE problem and plot temperature distribution
clc
close all

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

switch element_type
    case "LinearTriangle"
        nquadpoints = 1;
    case "BilinearQuadrilateral"
        nquadpoints = 2;
end


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
        error("conductivity cant be zero")
    end
    D = k * eye(2);
    
    % Solving Ke and fe
    switch element_type
        case "LinearTriangle"
            [Ke, fe] = linear_triangle_heat_element(D, 0, coords, nquadpoints);
        case "BilinearQuadrilateral"
            [Ke, fe] = bilinear_quadrilateral_heat_element(D, 0, coords, nquadpoints);
        otherwise
            error("bruh you dumb a shiieet")
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
neumann_dofs = node_sets{"hole"};

% Solve the equation system via partitioning
cdofs = constrained_dofs; % Dofs with known temperature, so constrain them!

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
            error("bruh you dumb a shiieet")
    end


end

% f(neumann_dofs) = q_2;
f = f + f_neumann;
fdofs = setdiff(1:ndofs, cdofs);    % The rest are unknown -> free to change

% Solve for the unknown dofs: af = Kff^{-1} * (ff - Kfc * ac)
af = K(fdofs, fdofs) \ (f(fdofs) - K(fdofs, cdofs) * a(cdofs));

a(fdofs) = af;

figure()
hold on
ae = get_element_dofs(element_nodes, a, 1);
[ex, ey] = get_element_coordinates(element_nodes, node_coordinates, element_type);
fill(ex, ey, ae);
axis equal;
colorbar

% ==========================================================
% calculate and plot heat flux in the center of each element
% ==========================================================

% Heat flux matrix
q = zeros(2, number_elements);

for i = 1:number_elements

    element_node_ids = element_nodes(:, i);          % Element node ids
    coords = node_coordinates(:, element_node_ids);  % Element node coordinates
    
    % Conductivity matrix
    k = k_el(i);
    if k == 0
        error("conductivity cant be zero")
    end
    D = k * eye(2);

    
    switch element_type
        case "LinearTriangle"
            q(:,i) = linear_triangle_heat_flux(D, coords, ae(:,i), 1);
        case "BilinearQuadrilateral"
            q(:,i) = bilinear_quadrilateral_heat_flux(D, coords, ae(:,i), 1);
        otherwise
            error("bruh you dumb a shiieet")
    end
end

[ex, ey] = get_element_coordinates(element_nodes, node_coordinates, element_type);

X = mean(ex);
Y = mean(ey);
U = q(1,:);
V = q(2,:);

quiver(X, Y, U, V, 'r', AutoScaleFactor=1.1)
hold off

%%
clc
% =============================================
% Total heat flux 
% =============================================

% a)
element_floortiles = facet_sets{"tilesurface"};
q_floortiles = q(:,element_floortiles(1,:));
facet_nr_floortiles = element_floortiles(2,:);

for i = 1:

switch element_type
    case "LinearTriangle"
        calculate_jacobian()
        q(:,i) = linear_triangle_heat_flux(D, coords, ae(:,i), 1);
        triangle_facet_weighted_normal()
    case "BilinearQuadrilateral"
        q(:,i) = bilinear_quadrilateral_heat_flux(D, coords, ae(:,i), 1);
    otherwise
        error("bruh you dumb a shiieet")
end


%%


%fprintf("Total flux into the room from the floor: %.1f W/m2\n", total_flux_into_room)
%fprintf("Total flux into the wall from the room: %.1f W/m2\n", total_flux_into_wall)
%fprintf("Total flux from the wall out into outside: %.1f W/m2\n", total_flux_out)