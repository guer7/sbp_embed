% Converts a gridfunction to a matrix
% Takes a grid function and a structured grid.
function F = funcToMatrix(g, gf)
    F = reshapeRowMaj(gf, g.size());
end