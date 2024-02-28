% Creates a curvilinear grid of dimension length(m).
% over the logical domain xi_lim, eta_lim, ...
% The grid point distribution is based on the boundary
% optimized SBP operators, passed as arguments, order
% and opt.
% Examples:
%   g = grid.boundaryOptimizedCurvilinear(mapping, [mx, my], xlim, ylim, order, opt)
%   g = grid.boundaryOptimizedCurvilinear(mapping, [10, 15], {0,1}, {0,2}, 4) - defaults to 'accurate' stencils
%   g = grid.boundaryOptimizedCurvilinear(mapping, [10, 15], {0,1}, {0,2}, 4, 'minimal')
function g = GLCurvilinear(mapping, m, varargin)
    n = length(m);
    % TODO: The input parameter check is also present in boundaryOptimizedGrid.m,
    % Consider refactoring it.
        
    % Check that parameters matches dimensions
    matchingParams = iscell([varargin{1:n}]) && isfloat([varargin{n+1}]);
    assert(matchingParams,'grid:GLCurvilinear:NonMatchingParameters','The number of parameters per dimensions do not match.');


    X = {};
    h = [];
    for i = 1:n
        try
            [X{i},h(i)] = sbp.grid.GLGrid(varargin{i},m(i),varargin{n+1});
        catch exception % Propagate any errors in the grid generation functions.
            msgText = getReport(exception);
            error('grid:boundaryOptimizedCurvilinear:InvalidParameter',msgText)
        end
    end

    g = grid.Curvilinear(mapping, X{:});
    g.logic.h = h;
end