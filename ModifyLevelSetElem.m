function [ffNew] = ModifyLevelSetElem(ffOld)
%This routine modifies the original level-set function such that the
%zero-level builds a straight line in the reference element(using 
%bi-linear shape functions for the interpolation of the level-set values
%from the nodes into the domain). This is the case if all four element
%nodes togeter with their level-set value build a plane. That is, 
%(r1, s1,f1), (r2, s2,f2),(r3, s3, f3), (r4, s4, f4) lie on a plane where
%fi are the level-set values.

% Important: This modification of the original level-set function does not
% effect the intersection points of the zero-level with the element edges.

% Important: This modification should be used with the sign-enrichment
% only, for the abs-enrichment it is not fully consistent.

% Definition of edge-number and corresponding nodes:
%    3    .3.    4 
%    o-----------o
%    |           |
%    |           |
% .4.|           |.2.
%    |           |   
%    |           |     .x. is the edge number.
%    o-----------o
%    1    .1.    2

ffNew = ffOld;
SignVect = sign(ffOld);
if isempty(find(SignVect == 0)) ==0
    error('Level set function is exatly zero at a node.')
end

if min(SignVect) == max(SignVect)
    return % Nothing to do if element is not cut
end

[CutEdge1, CutEdge2] = GetCutEdges(ffOld);

% Compute modified level-set values.
if CutEdge1 == 1 && CutEdge2 == 2
    ffNew(3) = ffNew(2)/ffOld(2) * ffOld(3);
    ffNew(4) = ffNew(1)-ffNew(2)+ffNew(3);
elseif CutEdge1 == 1 && CutEdge2 == 3
    ffNew(3) = ffOld(3) * (ffNew(1)-ffNew(2)) / (ffOld(4)-ffOld(3));
    ffNew(4) = ffOld(4) * (ffNew(1)-ffNew(2)) / (ffOld(4)-ffOld(3));
elseif CutEdge1 == 1 && CutEdge2 == 4
    ffNew(4) = ffNew(1)/ffOld(1) * ffOld(4);
    ffNew(3) = ffNew(4)-ffNew(1)+ffNew(2);
elseif CutEdge1 == 2 && CutEdge2 == 1
    ffNew(1) = ffNew(2)/ffOld(2) * ffOld(1);
    ffNew(4) = ffNew(1)-ffNew(2)+ffNew(3);
elseif CutEdge1 == 2 && CutEdge2 == 3
    ffNew(4) = ffNew(3)/ffOld(3) * ffOld(4);
    ffNew(1) = ffNew(2)-ffNew(3)+ffNew(4);
elseif CutEdge1 == 2 && CutEdge2 == 4
    ffNew(1) = ffOld(1) * (ffNew(2)-ffNew(3)) / (ffOld(1)-ffOld(4));
    ffNew(4) = ffOld(4) * (ffNew(2)-ffNew(3)) / (ffOld(1)-ffOld(4));
elseif CutEdge1 == 3 && CutEdge2 == 1
    ffNew(1) = ffOld(1) * (ffNew(4)-ffNew(3)) / (ffOld(1)-ffOld(2));
    ffNew(2) = ffOld(2) * (ffNew(4)-ffNew(3)) / (ffOld(1)-ffOld(2));
elseif CutEdge1 == 3 && CutEdge2 == 2
    ffNew(2) = ffNew(3)/ffOld(3) * ffOld(2);
    ffNew(1) = ffNew(2)-ffNew(3)+ffNew(4);
elseif CutEdge1 == 3 && CutEdge2 == 4
    ffNew(1) = ffNew(4)/ffOld(4) * ffOld(1);
    ffNew(2) = ffNew(3)-ffNew(4)+ffNew(1);
elseif CutEdge1 == 4 && CutEdge2 == 1
    ffNew(2) = ffNew(1)/ffOld(1) * ffOld(2);
    ffNew(3) = ffNew(4)-ffNew(1)+ffNew(2);
elseif CutEdge1 == 4 && CutEdge2 == 2
    ffNew(2) = ffOld(2) * (ffNew(1)-ffNew(4)) / (ffOld(2)-ffOld(3));
    ffNew(3) = ffOld(3) * (ffNew(1)-ffNew(4)) / (ffOld(2)-ffOld(3));
elseif CutEdge1 == 4 && CutEdge2 == 3
    ffNew(3) = ffNew(4)/ffOld(4) * ffOld(3);
    ffNew(2) = ffNew(3)-ffNew(4)+ffNew(1);
else
    error('Internal error.')
end

% Check the correctness of the result: Are all points (r_i, s_i, f_i) in the
% same plane and are the intersection points (ra, sa) and (rb, sb) unchanged?
CheckResult(ffOld, ffNew)

% % Plot situation.
% reset(cla), reset(clf), hold on
% patch([-1 1 1 -1], [-1 -1 1 1], 'y')
% a = patch([-1 1 1 -1], [-1 -1 1 1], ffOld, 'b');
% set(a, 'FaceAlpha', 0.2)
% a = patch([-1 1 1 -1], [-1 -1 1 1], ffNew, 'r');
% set(a, 'FaceAlpha', 0.2)
% a(1) = text(-1.1,-1.1, sprintf('%1.3f', ffOld(1)));
% a(2) = text( 1.1,-1.1, sprintf('%1.3f', ffOld(2)));
% a(3) = text( 1.1, 1.1, sprintf('%1.3f', ffOld(3)));
% a(4) = text(-1.1, 1.1, sprintf('%1.3f', ffOld(4)));
% set(a, 'FontSize', 20)
% view(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckResult(ffOld, ffNew)

    % Check that (r1, s1, fNew1), (r2, s2, fNew2), (r3, s3, fNew3), (r4, s4, fNew4) build
    % a plane.
    r1 = -1; s1 = -1;
    r2 =  1; s2 = -1;
    r3 =  1; s3 =  1;
    r4 = -1; s4 =  1;
    A = [r1-r2 r1-r3; s1-s2 s1-s3];
    b = [r4-r1; s4-s1];
    Sol = A\b;
    f4Check = ffNew(1) + Sol(1) * (ffNew(1)-ffNew(2)) + Sol(2)*(ffNew(1)-ffNew(3));
    if abs(ffNew(4)-f4Check) > 1.e-12
        error('Internal error: The new level set values build no plane!')
    end

    % Check that the intersection points of the disc. with the element edge are un-changed.
    if abs(sum(sign(ffOld)-sign(ffNew)))
        error('Internal error: The new level set function has different signs than the old one!')
    end

    Count = 0;
    if sign(ffOld(1)) ~= sign(ffOld(2))
        Count = Count + 1;
        rOld(Count) = Interpolate(r1, r2, ffOld(1), ffOld(2));
        sOld(Count) = Interpolate(s1, s2, ffOld(1), ffOld(2));
        rNew(Count) = Interpolate(r1, r2, ffNew(1), ffNew(2));
        sNew(Count) = Interpolate(s1, s2, ffNew(1), ffNew(2));
    end 
    if sign(ffOld(2)) ~= sign(ffOld(3))
        Count = Count + 1;
        rOld(Count) = Interpolate(r2, r3, ffOld(2), ffOld(3));
        sOld(Count) = Interpolate(s2, s3, ffOld(2), ffOld(3));
        rNew(Count) = Interpolate(r2, r3, ffNew(2), ffNew(3));
        sNew(Count) = Interpolate(s2, s3, ffNew(2), ffNew(3));
    end 
    if sign(ffOld(3)) ~= sign(ffOld(4))
        Count = Count + 1;
        rOld(Count) = Interpolate(r3, r4, ffOld(3), ffOld(4));
        sOld(Count) = Interpolate(s3, s4, ffOld(3), ffOld(4));
        rNew(Count) = Interpolate(r3, r4, ffNew(3), ffNew(4));
        sNew(Count) = Interpolate(s3, s4, ffNew(3), ffNew(4));
    end 
    if sign(ffOld(4)) ~= sign(ffOld(1))
        Count = Count + 1;
        rOld(Count) = Interpolate(r4, r1, ffOld(4), ffOld(1));
        sOld(Count) = Interpolate(s4, s1, ffOld(4), ffOld(1));
        rNew(Count) = Interpolate(r4, r1, ffNew(4), ffNew(1));
        sNew(Count) = Interpolate(s4, s1, ffNew(4), ffNew(1));
    end 

    if sum(sum(abs([rOld; sOld] - [rNew; sNew]))) > 1.e-12
        error('Internal error: The intersection points have changed between old and new level set!')
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function[CutEdge1, CutEdge2] = GetCutEdges(ffElem)
        
        SignVect = sign(ffElem) ;
        if isempty(find(SignVect == 0)) == 0
            error ('Level-set is exactly zero at a node!')
        end
        
        if min(SignVect) == max(SignVect)
            error('This is not a cut element')
        end
        
        Count = 0;
        if SignVect(1) ~= SignVect(2)
            Count = Count +1;
            CutEdge(Count, 1) = 1;
        end
         if SignVect(2) ~= SignVect(3) % Disc. is between node 2 and 3.
           Count = Count + 1;
           CutEdge(Count, 1) = 2;
         end 
          if SignVect(3) ~= SignVect(4) % Disc. is between node 3 and 4.
        Count = Count + 1;
        CutEdge(Count, 1) = 3;
        end
             if SignVect(4) ~= SignVect(1) % Disc. is between node 4 and 1.
        Count = Count + 1;
        CutEdge(Count, 1) = 4;
    end

    if Count == 4
        ffElem
        error('This case for the level-set function is not allowed (disc. not uniquely defined!')
    elseif Count ~= 2
        ffElem
        error('Internal error!')
    end
    
    CutEdge1 = CutEdge(1);
    CutEdge2 = CutEdge(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xStar] = Interpolate(x1, x2, f1, f2)

    xStar = x1 + (x2-x1) * f1 / (f1-f2);
    
end