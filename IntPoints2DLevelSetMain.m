function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSetMain(xxElem, ...
    yyElem, ffElem, NodesAct1, NodesAct2, CrackTip, nQ)

% Get the Gauss-Legendre integration points. Note that almost-polar
% integration is used in the crack-tip element.

if length(find(NodesAct1 == 0)) ==4 && length(find(NodesAct2==0))==4
    % Standard finite element
    % patch(xxElem,yyElem, 'b')       Plot the standard finite element
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DRefElemQuad(2);
    % xxIntRef  and yyIntRef  Guass-Legendre Integration Point
    % wwIntRef     Guass-Legendre Integration Point   
elseif isempty(find(xxElem>CrackTip.xx))==0 && isempty(find(xxElem<CrackTip.xx)) == 0 &&...
        isempty(find(yyElem>CrackTip.yy))==0 && isempty(find(yyElem<CrackTip.yy)==0)
    % Crack-Tip element.
    % patch(xxElem, yyElem, 'r')
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem, yyElem, ffElem, 0 , 0, 5);
    %[xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
else %All other elements
    [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DLevelSet(ffElem, nQ);
end
