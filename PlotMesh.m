function PlotMesh(Mesh, xx, yy, ElemNum)

% Plot the mesh.

reset(cla), reset(clf), hold on

for  CurrElem = 1 : ElemNum
     Nodes = Mesh(CurrElem, :);
     patch(xx(Nodes), yy(Nodes), -0.005*[0.5 0.5 0.5 0.5], 'w')
end

axis equal
axis([min(xx)-0.1   max(xx)+0.1     min(yy)-0.1     max(yy)+0.1])