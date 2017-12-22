function [N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(...
    xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQ)

% Compute standard bi-linear shape functions at integration points in the
% real element.

% Initialization.

N = zeros(4, nQ);
dNdx = zeros(4, nQ);
dNdy = zeros(4, nQ);

xxInt = zeros(1, nQ);
yyInt = zeros(1,nQ);
wwInt = zeros(1, nQ);

% Define N, dNdr, dNds at the integration points in the REFERENCE element.
