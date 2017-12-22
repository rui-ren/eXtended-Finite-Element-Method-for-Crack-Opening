function [NodeNum, ElemNum, xx, yy, Mesh] = GetMesh (lx, ly, nx, ny)
% Create a quadrilateral mesh in a domain of size[0, lx] * [0, ly]
% nx elements in x-direction and ny element in y-direction.
% Set number to every element from bottom to top and left to right
% y
% ^
% |
% |
% + ----------> 
%            x

%%%%%%%%%%%%%%%%Geometry of the Domain%%%%%%%%%%%%%%
HelpVec = linspace(ly, 0, ny+1);           % 
NodeNum = (nx + 1)*(ny + 1);               % the whole Node in the domain
ElemNum = nx * ny;                         % the whole Element in the domain

% Set node   Initialization the node
xx = zeros(NodeNum, 1);
yy = zeros(NodeNum, 1);

for i = 1 : nx +1          % loop the node on the x-coordination
    xx((i - 1)*(ny + 1)+1 : i *(ny +1)) = (i - 1)*(lx / nx);          % Node x-value--- iteration the every line
    yy((i - 1)*(ny + 1)+1 : i *(ny + 1)) = HelpVec;                   % Node y-value----iteration the every line
end


%%%%%%%%%%%%%%%%%%Mesh of the Domain %%%%%%%%%%%%%%%%

%Define element
Mesh = zeros(ElemNum, 4);
for i = 1: nx
    for j = 1 : ny
        CurrElem = (i -1)*ny + j;
        c = (i -1)*(ny + 1) + j;                                                                  % the top - left Node for every Element
        Mesh(CurrElem, :) = [c+1    c+1+(ny+1)    c+ (ny+1)     c];      % assemble the Node to every Element
    end
end

%%%%%%%%%%%%%%%%%%Perturb inner nodes %%%%%%%%%%%%%%%%

%Add the Perturbation Theory to the mesh generation method.
CasePerturb = 0;
if CasePerturb == 1
    ScalPerturb = 0.2;   % Scales the magnitude of the perturbation.
    xMaxPerturb = 1/2*(lx/nx);
    yMaxPerturb = 1/2*(ly/ny);
    rand('state', 0)
    xPerturb = ScalPerturb * xMaxPerturb * (2*rand(NodeNum, 1) -1);
    yPerturb = ScalPerturb * yMaxPerturb * (2*rand(NodeNum,1) -1);
    Bound = find(xx==0|xx==lx|yy==0|yy==ly);      % find the boundary point in the domain
    PosInner = setdiff(1:NodeNum, Bound);                  % eliminate the boundary point in the system
    
    xx(PosInner) = xx(PosInner) + xPerturb(PosInner);
    yy(PosInner) = yy(PosInner) + yPerturb(PosInner);   % perturbation theory.
end