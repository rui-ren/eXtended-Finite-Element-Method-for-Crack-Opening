function [xxIntRef, yyIntRef, wwIntRef] = IntPoints2DCrackTipElem(xxElem,...
    yyElem, ffElem, xxTip, yyTip, nQ1D)

% Divide the element containing the crack tip into six integration point
% each triangle contains the crack tip as a node.

% Input data are the real element coordinates, the level set function and
% the crack tip position. Note, that the level set is a straight line in
% the reference element. But the crack tip position is defined in the real
% domain

% Extract the cut element segments of this element.
SignVect = sign(ffElem);

% ??????????
% ???????? Guass Integration Point
Count = 0;
if SignVect(1) ~= SignVect(2)
    Count = Count + 1;
    CutSegm(Count) = 1;
    xxS(Count) = Interpolate(-1,  1, ffElem(1), ffElem(2));
    yyS(Count) = Interpolate(-1, -1, ffElem(1), ffElem(2));
end
if SignVect(2) ~= SignVect(3)
    Count = Count +1;
    CutSegm(Count) = 2;
    xxS(Count) = Interpolate(  1, 1,  ffElem(2), ffElem(3));
    yyS(Count) = Interpolate(-1, 1, ffElem(3), ffElem(3));
end
if SignVect(3) ~= SignVect(4)
    Count = Count +1;
    CutSegm(Count) = 3;
    xxS(Count) = Interpolate(1, -1, ffElem(3), ffElem(4));
    yyS(Count) = Interpolate(1,  1, ffElem(3), ffElem(4));
end
if SignVect(4) ~= SignVect(1)
    Count = Count +1;
    CutSegm(Count) = 4;
    xxS(Count) = Interpolate(-1, -1,ffElem(4), ffElem(1));
    yyS(Count) = Interpolate(1, -1, ffElem(4), ffElem(1));
    
    if Count ~=2
        error('Internal Error.')
    end
    
    % Construct the six integration triangles...
    % Number the nodes.
    
    Count = 1;
    for i = 1 : 4
        HelpNodes(i + Count - 1) = i;
        if Count <= 2
            if CutSegm(Count) == i
                HelpNodes(i + Count) = 4+Count;
                Count = Count + 1;
            end
        end
    end
    
    % ... Set Node Coordinates.
    
    [xxTipRef, yyTipRef, CaseIntPointInElem] = ...
        ProjectRealToRefElem(xxElem, yyElem, xxTip, yyTip);
    
    if CaseIntPointInElem == 0
        error('Crack tip is not in supposed crack-tip element!')
    end
    
    xxTri = [[-1 1 1 -1], xxS(1), xxS(2), yyTipRef];
    yyTri = [[-1 -1 1 1], yyS(1), yyS(2), yyTipRef ];
    
    for i = 1 : 6
        MeshTri(i, 1) = 7;
        MeshTri(i, 2) = HelpNodes(i);
            if i == 6
                MeshTri(i, 3) = HelpNodes(1);
            else
                MeshTri(i, 3) = HelpNodes(i +1);
            end
    end
    
     % Set integration points in each triangle according to almost polar integration.
    xxIntRef = zeros(1, 6*nQ1D*nQ1D);
    yyIntRef = zeros(1, 6*nQ1D*nQ1D);
    wwIntRef = zeros(1, 6*nQ1D*nQ1D);
    [xxIntRefTri, yyIntRefTri, wwIntRefTri] = IntPoints2DRefElemTri(nQ1D);
    for i = 1 : 6
    NodesTri = MeshTri(i, :);
    xxElemTri = xxTri(NodesTri);
    yyElemTri = yyTri(NodesTri);
    [xxInt, yyInt, wwInt] = IntPoints2DRealElemTri(xxElemTri, yyElemTri, ...
        xxIntRefTri, yyIntRefTri, wwIntRefTri, nQ1D*nQ1D);
    nQ2D = nQ1D * nQ1D;
    xxIntRef((i-1)*nQ2D+1 : i*nQ2D) = xxInt;
    yyIntRef((i-1)*nQ2D+1 : i*nQ2D) = yyInt;
    wwIntRef((i-1)*nQ2D+1 : i*nQ2D) = wwInt;
    end
   
end

function [xStar] = Interpolate(x1, x2, f1, f2)
        xStar = x1 + (x2 - x1)*f1/(f1 - f2);
    end