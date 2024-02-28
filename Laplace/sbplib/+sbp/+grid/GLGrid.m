% Computes the grid points x and grid spacing h used by the boundary optimized SBP operators
% with minimal number of non-equidistant boundary points, presented in
% 'Boundary optimized diagonal-norm SBP operators - Mattsson, Almquist, van der Weide 2018'.
%
% lim - cell array with domain limits
% N - Number of grid points
% order - order of accuracy of sbp operator.
function [x,h] = GLGrid(lim,N,order)
assert(iscell(lim) && numel(lim) == 2,'The limit should be a cell array with 2 elements.');
L = lim{2} - lim{1};
assert(L>0,'Limits must be given in increasing order.');

switch order
    case 3
        assert(N == 4, 'The number of grid points must be order + 1')
        d1 = 1 - sqrt(5)/5;
        x = 0.5*[-1 (d1-1) 1-d1 1]' + 0.5;
    case 4
        assert(N == 5, 'The number of grid points must be order + 1')
        d1 = 1 - sqrt(21)/7;
        x = 0.5*[-1 (d1-1) 0 1-d1 1]' + 0.5;
    case 5
        assert(N == 6, 'The number of grid points must be order + 1')
        d1 = 1 - sqrt(147 + 42*sqrt(7))/21;
        d2 = sqrt(147 + 42*sqrt(7))/21 - sqrt(147 - 42*sqrt(7))/21;
        x = 0.5*[-1 (d1-1) (d1+d2-1) 1-(d1+d2) 1-d1 1]' + 0.5;
    case 6
        assert(N == 7, 'The number of grid points must be order + 1')
        d1 = -sqrt(990 - 66*sqrt(165))/66 + 1 - sqrt(990 + 66*sqrt(165))/66;
        d2=sqrt(990 - 66*sqrt(165))/33;
        x = 0.5*[-1 (d1-1) (d1+d2-1) 0 1-(d1+d2) 1-d1 1]' + 0.5;
    case 7
        assert(N == 8, 'The number of grid points must be order + 1')
        x = 0.5*[-0.10e1; -0.87174014850960661533e0; -0.59170018143314230214e0; -0.20929921790247886873e0; 0.20929921790247886873e0; 0.59170018143314230214e0; 0.87174014850960661533e0; 0.10e1;] + 0.5;
    case 9
        assert(N == 10, 'The number of grid points must be order + 1')
        x = 0.5*[-1; -0.919533908166458813828932660822e0; -0.738773865105505075003106174860e0; -0.477924949810444495661175092731e0; -0.165278957666387024626219765958e0; 0.165278957666387024626219765958e0; 0.477924949810444495661175092731e0; 0.738773865105505075003106174860e0; 0.919533908166458813828932660822e0; 1;] + 0.5;
end

x = lim{1} + (lim{2} - lim{1})*x;
h = L/(N-1);
end
