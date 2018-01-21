function [ElemNum1D, BoundaryElems1D, CorrBoundaryElems2D] = ...
    DetectBoundaryElems(NodeNum, ElemNum, xx, yy, Mesh)

% Construct 1D boundary element data. Therefore, the element edges on the
% boundary are extracted. ElemNum 1D is the number of 1D elements, Boundary
% Elems1D stores, the element nodes corresponding to the global 2D mesh,
% and CorrBoundaryElems2D stores which 2D element belongs to each boundary
% element.

Bound = find (xx == min(xx) | xx == max(xx) | yy == min(yy) | yy == max(yy));

% Detect 2D boundary elements.

BoundaryElems2D = zeros(ElemNum, 1);

for m = 1 : ElemNum
    Nodes = Mesh(m, :);
    RemainingNodes = setdiff(Nodes, Bound);
    
   if length(RemainingNodes) < 4          % ???????
      BoundaryElems2D(m,1) = 1;          %  ????????
   end
end

    BoundaryElems2D = find(BoundaryElems2D==1);

    if isempty(BoundaryElems2D) == 1
        error('Internal error')
    end

% Extract 1D boundary elements.
    Segment1D = zeros(4,2);
    count = 0;

for m = 1 : length(BoundaryElems2D)
    
    Nodes = Mesh(BoundaryElems2D(m), :);
    
    Segment1D(1, :) = [Nodes(1)   Nodes(2)];
    Segment1D(2, :) = [Nodes(2)   Nodes(3)];
    Segment1D(3, :) = [Nodes(3)   Nodes(4)];
    Segment1D(4, :) = [Nodes(4)   Nodes(1)];
    
    for i = 1 : 4
        ElementBoundaryNodes = intersect(Segment1D(i, :), Bound);
        if length(ElementBoundaryNodes) == 2
            % It may still be possible, that the two nodes belong to different boundary parts.
            % Then, this is in fact an INNER element edge. Therefore, make the following check:
            if (xx(Segment1D(i, 1)) == min(xx) & xx(Segment1D(i, 2)) == min(xx)) | ...
               (xx(Segment1D(i, 1)) == max(xx) & xx(Segment1D(i, 2)) == max(xx)) | ...
               (yy(Segment1D(i, 1)) == min(yy) & yy(Segment1D(i, 2)) == min(yy)) | ...
               (yy(Segment1D(i, 1)) == max(yy) & yy(Segment1D(i, 2)) == max(yy))
                count = count + 1;
                BoundaryElems1D(count, :) = Segment1D(i, :);
                CorrBoundaryElems2D(count, 1) = BoundaryElems2D(m);
            end
        end
    end
end

ElemNum1D = count;
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

      
    
    
