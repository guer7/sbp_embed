% Setup closure and penalty matrices for several boundary conditions at once.
% Each bc is a struct with the fields
%  * type     -- Type of boundary condition
%  * boundary -- Boundary identifier
%  * data     -- A function_handle for a function which provides boundary data.(see below)
% Also takes S_sign which modifies the sign of the penalty function, [-1,1]
% Returns a closure matrix and a penalty matrices for each boundary condition.
%
% The boundary data function can either be a function of time or a function of time and space coordinates.
% In the case where it only depends on time it should return the data as grid function for the boundary.
% In the case where it also takes space coordinates the number of space coordinates should match the number of dimensions of the problem domain.
% For example in the 2D case: f(t,x,y).
function [closure, penalties, L] = closureSetup(diffOp, bcs, is_derv)
scheme.bc.verifyFormat(bcs, diffOp);

% Setup storage arrays
closure = {spzeros(size(diffOp)),spzeros(size(diffOp))};
penalties = cell(1, length(bcs));

% Collect closures and penalties
L = [];
for i = 1:length(bcs)
    switch bcs{i}.method
        case 'sat'
            if is_derv
                [localClosure, localpenalties] = diffOp.boundary_condition_dervs(bcs{i}.boundary, bcs{i}.type);
            else
                [localClosure, localpenalties] = diffOp.boundary_condition(bcs{i}.boundary, bcs{i}.type);
            end
            penalties{i} = localpenalties;
%             if ~isempty(diffOp.Pic)
%                localClosure = diffOp.Pic*localClosure*diffOp.Pic;
%                penalties{i} = diffOp.Pic*penalties{i};
%             end
            closure{1} = closure{1} + localClosure{1};
            closure{2} = closure{2} + localClosure{2};
        case 'proj'
            L = [L;
                diffOp.getBoundaryOperator(bcs{i}.type,bcs{i}.boundary)'];
    end
end
if ~isempty(L)
    L = remove_deps(L,diffOp.H);
end
end
