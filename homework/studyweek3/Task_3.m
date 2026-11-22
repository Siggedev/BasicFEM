%% LinearLine
clc
close all

fprintf("========================================================\n")
fprintf("======== LinearLine ========\n")
fprintf("========================================================\n")
[element_nodes, node_coordinates, facet_sets, node_sets] = generate_mesh( ...
    'LinearLine', 2, 0.5, 2)

fprintf("======== facet_sets ========")
facet_left = facet_sets{'left'}
facet_left = facet_sets{'right'}

%% QuadraticLine
clc
close all
fprintf("========================================================\n")
fprintf("======== QuadraticLine ========\n")
fprintf("========================================================\n")

[element_nodes, node_coordinates, facet_sets, node_sets] = generate_mesh( ...
    'QuadraticLine', 2, 0.5, 2)

fprintf("======== facet_sets ========")
facet_left = facet_sets{'left'}
facet_left = facet_sets{'right'}
%% LinearTriangle
clc
close all
fprintf("========================================================\n")
fprintf("======== LinearTriangle ========\n")
fprintf("========================================================\n")

[element_nodes, node_coordinates, facet_sets, node_sets] = generate_mesh( ...
    'LinearTriangle', [5; 5])

[ex, ey] = get_element_coordinates(element_nodes, node_coordinates, 'LinearTriangle');
fill(ex, ey, 'w')
axis equal;
plot_fe_mesh(element_nodes, node_coordinates, "LinearTriangle")

%% Bilinear Quadrilateral
clc
close all
[element_nodes, node_coordinates, facet_sets, node_sets] = generate_mesh( ...
    'BilinearQuadrilateral', [5,5])
[ex, ey] = get_element_coordinates(element_nodes, node_coordinates, 'BilinearQuadrilateral');
fill(ex, ey, 'w')
axis equal;
plot_fe_mesh(element_nodes, node_coordinates, "BilinearQuadrilateral")

%% Read meshes
clc
close all

hw3_testmesh_LinearTriangle_path = fullfile("homework", "studyweek3/hw3_testmesh_LinearTriangle.inp");

[element_nodes, node_coordinates, element_sets, facet_sets, node_sets, element_type]= read_mesh(hw3_testmesh_LinearTriangle_path);
[ex, ey] = get_element_coordinates(element_nodes, node_coordinates, 'LinearTriangle');

% axis equal;
% fill(ex, ey, 'w') 
% fig = plot_fe_mesh(element_nodes, node_coordinates, 'LinearTriangle');

el_data = zeros(size(element_nodes, 2), 1);
el_data(element_sets{'floorpanel'}) = 1;
fill(ex, ey, el_data)
axis equal;
colorbar();
