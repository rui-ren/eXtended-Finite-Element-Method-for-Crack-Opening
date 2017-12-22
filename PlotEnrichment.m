function PlotEnrichment(Mesh, xx, yy, ElemsEnriched, NodesEnriched)

%Plot the enriched elements and nodes

for i = 1: length(ElemsEnriched)
    CurrElem = ElemsEnriched(i);
    Nodes = Mesh(CurrElem,:);
    patch(xx(Nodes), yy(Nodes), -0.004*[1 1 1 1], 'r')
end

plot3(xx(NodesEnriched), yy(NodesEnriched), -0.001*ones(length(NodesEnriched), 1), 'ko')

end

