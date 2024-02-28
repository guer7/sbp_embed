% Numbering of boundaries on reference domain.
function idx = bound2idx(bound)
switch bound
    case 'w'
        idx = 1;
    case 'e'
        idx = 2;
    case 's'
        idx = 3;
    case 'n'
        idx = 4;
end