% Creates a curvilinear grid of dimension length(m).
% over the logical domain xi_lim, eta_lim, ...
% The grid point distribution is based on the boundary
% optimized SBP operators, passed as arguments, order
% and opt.
% Examples:
%   g = grid.boundaryOptimizedCurvilinear(mapping, [mx, my], xlim, ylim, order, opt)
%   g = grid.boundaryOptimizedCurvilinear(mapping, [10, 15], {0,1}, {0,2}, 4) - defaults to 'accurate' stencils
%   g = grid.boundaryOptimizedCurvilinear(mapping, [10, 15], {0,1}, {0,2}, 4, 'minimal')
function g = boundaryOptimizedCurvilinear(mapping, m, varargin)
    n = length(m);
        
    % Check that parameters matches dimensions
    matchingParams = false;
    if length(varargin) == n+1 % Minimal number of arguments
            matchingParams = iscell([varargin{1:n}]) && ...
                             isfloat([varargin{n+1}]);
    elseif length(varargin) == n+2 % Stencil options supplied
            matchingParams = iscell([varargin{1:n}]) && ...
                             isfloat([varargin{n+1}]) && ...
                             ischar([varargin{n+2}]);
    end
    assert(matchingParams,'grid:boundaryOptimizedCurvilinear:NonMatchingParameters','The number of parameters per dimensions do not match.');

    % Check that stencil options are passed correctly (if supplied)
    if length(varargin) == n+2 % Stencil options supplied
        availabe_opts = ["Accurate","accurate","A","acc","Minimal","minimal","M","min"];
        assert(any(varargin{n+2} == availabe_opts), ...
            'grid:boundaryOptimizedCurvilinear:InvalidOption',"The operator option must be 'accurate' or 'minimal.'");
    else %If not passed, populate varargin with default option 'accurate'
        varargin(n+2) = {'accurate'};
    end

    % Specify generating function
    switch varargin{n+2}
        case {'Accurate','accurate','A','acc'}
            gridgenerator = @sbp.grid.accurateBoundaryOptimizedGrid;
        case {'Minimal','minimal','M','min'}
            gridgenerator = @sbp.grid.minimalBoundaryOptimizedGrid;
    end

    X = {};
    h = [];
    for i = 1:n
        try
            [X{i},h(i)] = gridgenerator(varargin{i},m(i),varargin{n+1});
        catch exception % Propagate any errors in the grid generation functions.
            msgText = getReport(exception);
            error('grid:boundaryOptimizedCurvilinear:InvalidParameter',msgText)
        end
    end

    g = grid.Curvilinear(mapping, X{:});
    g.logic.h = h;
end