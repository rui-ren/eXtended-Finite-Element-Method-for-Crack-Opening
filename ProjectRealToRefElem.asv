function [r, s, CaseIntPointInElem] = ProjectRealToRefElem(xElem, yElem, xInt, yInt)

% Project Current Integration point into[-1, 1] * [-1 1] reference element.
% Note: This routine does not consider all special cases where the
% transformation between the real and reference element is singular with
% the general procedure.

x1 = xElem(1);  x2 = xElem(2);   x3 = xElem(3);  x4 = xElem(4);
y1 = yElem(1);  y2 = yElem(2);   y3 = yElem(4);  y4 = yElem(4);

EPS = 1.e-10;
if x1==x2 && x3==x4
    s = 2 * (xInt-x1) / (x3-x1) - 1;
    r = (4*yInt - ( (1-s)*y1+(1-s)*y2+(1+s)*y3+(1+s)*y4 ) ) / ...
                  (-(1-s)*y1+(1-s)*y2+(1+s)*y3-(1+s)*y4 );
    if abs(s)>=1+EPS || abs(r)>=1+EPS %Integration point not in element!
        CaseIntPointInElem = 0;
        return
    end   
elseif x2==x3 && x4==x1
    r = 2 * (xInt-x1) / (x2-x1) - 1;
    s = (4*yInt - ( (1-r)*y1+(1+r)*y2+(1+r)*y3+(1-r)*y4 ) ) / ...
                  (-(1-r)*y1-(1+r)*y2+(1+r)*y3+(1-r)*y4 );
    if abs(s)>=1+EPS || abs(r)>=1+EPS %Integration point not in element!
        CaseIntPointInElem = 0;
        return
    end   
elseif y1==y2 && y3==y4
    s = 2 * (yInt-y1) / (y3-y1) - 1;
    r = (4*xInt - ( (1-s)*x1+(1-s)*x2+(1+s)*x3+(1+s)*x4 ) ) / ...
                  (-(1-s)*x1+(1-s)*x2+(1+s)*x3-(1+s)*x4 );
    if abs(s)>=1+EPS || abs(r)>=1+EPS %Integration point not in element!
        CaseIntPointInElem = 0;
        return
    end   
elseif y2==y3 && y4==y1
    r = 2 * (yInt-y1) / (y2-y1) - 1;
    s = (4*xInt - ( (1-r)*x1+(1+r)*x2+(1+r)*x3+(1-r)*x4 ) ) / ...
                  (-(1-r)*x1-(1+r)*x2+(1+r)*x3+(1-r)*x4 );
    if abs(s)>=1+EPS || abs(r)>=1+EPS %Integration point not in element!
        CaseIntPointInElem = 0;
        return
    end   
else
    
    % Transform all others like this.
    px = x1 + x2 + x3 + x4;
    ux =-x1 + x2 + x3 - x4;
    vx = x1 - x2 + x3 - x4;
    tx =-x1 - x2 + x3 + x4;

    py = y1 + y2 + y3 + y4;
    uy =-y1 + y2 + y3 - y4;
    vy = y1 - y2 + y3 - y4;
    ty =-y1 - y2 + y3 + y4;
    
    help1=-4*xInt*uy+px*uy+4*yInt*ux-py*ux;
    help2=-4*xInt*vy+px*vy+tx*uy+4*yInt*vx-py*vx-ty*ux;
    help3=tx*vy-ty*vx;

    p = help2/help3;
    q = help1/help3;

    s1 = -1/2*p+1/2*(p^2-4*q)^(1/2);
    r1 = (4*xInt-px-s1*tx)/(ux+s1*vx);

    s2 = -1/2*p-1/2*(p^2-4*q)^(1/2);
    r2 = (4*xInt-px-s2*tx)/(ux+s2*vx);

    if abs(s1)<=1+EPS && abs(r1)<=1+EPS
        r=r1;
        s=s1;
    elseif abs(s2)<=1+EPS && abs(r2)<=1+EPS
        r=r2;
        s=s2;
    else % Integration point not in element!
        CaseIntPointInElem = 0;
        return
    end   

end

CaseIntPointInElem=1;
    
    
    
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
