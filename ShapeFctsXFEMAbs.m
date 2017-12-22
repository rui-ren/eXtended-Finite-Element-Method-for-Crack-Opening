
function [N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt] = ...
    ShapeFctsXFEMAbs(xxElem, yyElem, ffElem, NodesAct, xxIntRef, ...
    yyIntRef, wwIntRef, nQ)

% Shape Function and Enriched functions at integration points in the real
% element.

% Get standard FE shape function and level-set values at int. points.

[N, dNdx, dNdy, xxInt, yyInt, wwInt] = ShapeFctsStrdFEM(xxElem, yyElem, ...
    xxIntRef, yyIntRef, wwIntRef, nQ);

ffInt = N' * ffElem;

% Define the enriched functions, M, dMdx, dMdy for the enriched element
% nodes.

M = zeros(4, nQ);
dMdx = zeros(4, nQ);
dMdy = zeros(4, nQ);

if sum(NodesAct) == 0
    return                  % Go back if enriched is not active in this element.
end

for i = 1 : nQ
    
    % Strd. FEM shape fcts. at int.point.
    NInt = N( : , i);
    NxInt = dNdx(:, i);
    NyInt = dNdy(:, i);
    
    % Abs(Levelset) at int. point. and nodes.
    PsiElem = abs(ffElem);
    PsiInt = abs(ffInt(i));
    dPsidxInt = sign(ffInt(i)) * NxInt' * ffElem;
    dPsidyInt = sign(ffInt(i)) * NyInt' * ffElem;
    
    % Enrichment function at int. point.
    M(:, i)    = NInt.*(PsiInt-PsiElem);                % Strong discontinuity
    dMdx(:, i) = NxInt.*(PsiInt-PsiElem) + NInt.*dPsidxInt;
    dMdy(:, i) = NyInt.*(PsiInt-PsiElem) + NInt.*dPsidyInt;
    
end

Pos = find(NodesAct == 0);
M(Pos, :)    = 0;
dMdx(Pos, :) = 0;
dMdy(Pos, :) = 0;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    