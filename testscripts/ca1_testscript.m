% This test script should be used to validate your BasicFEM 
% implementations that are used in CA1. It will print "Pass" for all test
% cases that works correctly, and give more output if something fails. 
% 
% This is not written with the intention that you should modify this file, 
% but simply to be able to quickly validate if your implementations are
% correct, and if not, which functions that may be wrong. However, it might
% still be useful with this brief explanation of the datastructure:
% 
% The file "CA1.mat" contains a single variable, `data`, which is a 
% cell array of testcases. Each test case, contains two parts:
% 1) input: A dictionary describing the chosen input to generate the output
%        to be tested against. The fields are
%   * coords: Element coordinates
%   * D: Element conductivity matrix (Fourier's law)
%   * b: Element heat source
%   * ae: Element solution dofs
%   * qn: Boundary normal heat flux
%   * nqp: Number of quadrature points for the element
%   * nqp_nbc: Number of quadrature points for the neumann BCs.
%   * element_type: The element type used to generate data ("LinearTriangle"
%     or "BilinearQuadrilateral"
% 2) testdata: A dictionary describing the resulting output for different
%    variables with the given inputs. The fields are
%   * Ke: Element stiffness matrix
%   * fe: Element load vector 
%   * q: Calculated flux inside each quadrature point
%   * fe_nbc: Calculated Neumann contribution from each facet (cell array)
clc; clear all; close all;

testcases = load("CA1.mat").data;

for i = 1:length(testcases)
    testcase = testcases{i};
    input = testcase{1};
    testdata = testcase{2};
    coords = input{"coords"};
    D = input{"D"};
    b = input{"b"};
    ae = input{"ae"};
    qn = input{"qn"};
    nqp = input{"nqp"};
    nqp_nbc = input{"nqp_nbc"};
    
    num_failed = 0;
    switch input{"element_type"}
        case "BilinearQuadrilateral"
            [Ke, fe] = bilinear_quadrilateral_heat_element(D, b, coords, nqp);
            failed = test_values(input, testdata{"Ke"}, Ke, "Ke", i, "bilinear_quadrilateral_heat_element", "");
            num_failed = num_failed + failed;
            failed = test_values(input, testdata{"fe"}, fe, "fe", i, "bilinear_quadrilateral_heat_element", "");
            num_failed = num_failed + failed;
            q = bilinear_quadrilateral_heat_flux(D, coords, ae, nqp);
            failed = test_values(input, testdata{"q"}, q, "q", i, "bilinear_quadrilateral_heat_flux", "");
            num_failed = num_failed + failed;
            fe_nbc_truevals = testdata{"fe_nbc"};
            for facet_nr = 1:4
                fe_nbc = bilinear_quadrilateral_heat_neumann(qn, coords, facet_nr, nqp_nbc);
                fe_nbc_true = fe_nbc_truevals{facet_nr};
                failed = test_values(input, fe_nbc_true, fe_nbc, "fe_nbc", i, "bilinear_quadrilateral_heat_neumann", sprintf("Facet %u", facet_nr));
                num_failed = num_failed + failed;
            end
        case "LinearTriangle"
            [Ke, fe] = linear_triangle_heat_element(D, b, coords, nqp);
            failed = test_values(input, testdata{"Ke"}, Ke, "Ke", i, "linear_triangle_heat_element", "");
            num_failed = num_failed + failed;
            failed = test_values(input, testdata{"fe"}, fe, "fe", i, "linear_triangle_heat_element", "");
            num_failed = num_failed + failed;
            q = linear_triangle_heat_flux(D, coords, ae, nqp);
            failed = test_values(input, testdata{"q"}, q, "q", i, "linear_triangle_heat_flux", "");
            num_failed = num_failed + failed;
            fe_nbc_truevals = testdata{"fe_nbc"};
            for facet_nr = 1:3
                fe_nbc = linear_triangle_heat_neumann(qn, coords, facet_nr, nqp_nbc);
                fe_nbc_true = fe_nbc_truevals{facet_nr};
                failed = test_values(input, fe_nbc_true, fe_nbc, "fe_nbc", i, "linear_triangle_heat_neumann", sprintf("Facet %u", facet_nr));
                num_failed = num_failed + failed;
            end
        otherwise
            error("Something is wrong with the inputdata")
    end
    if num_failed == 0
        resultstring = "Pass";
    else
        resultstring = "Failed, please fix";
    end
    fprintf("Test case %u (%s, nqp = %u, nqp_nbc = %u): %s\n", ...
        i, input{"element_type"}, input{"nqp"}, input{"nqp_nbc"}, resultstring)
end


function failed = test_values(input, true_value, value, key, casenr, routine, extra_msg)
    if ~my_isapprox(value, true_value, 1e-8)
        eltype = input{"element_type"};
        fprintf("Test case %u failed for %s\n", casenr, key);
        fprintf("Wrong output from `%s`, \n", routine)
        fprintf("%s: nqp = %u, nqp_nbc = %u", eltype, input{"nqp"}, input{"nqp_nbc"})
        if ~isempty(extra_msg)
            fprintf("%s\n", extra_msg)
        end
        fprintf("Displaying your value and the correct value\n")
        disp(value)
        disp(true_value)
        failed = true;
    else
        failed = false;
    end
end
