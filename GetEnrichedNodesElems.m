function [ElemsEnriched, NodesEnriched] = GetEnrichedNodesElems(Mesh, ff)
% Get cut elements and corresponding nodes;
ElemsEnriched = find(min(sign(ff(Mesh))') ~= max(sign(ff(Mesh))'));  % enrich the element
ElemsEnriched = ElemsEnriched';                 %  transpose the element.
NodesEnriched = Mesh(ElemsEnriched, :);    % enrich the element in the cracks, the Node need be enriched
NodesEnriched = NodesEnriched(:);
NodesEnriched = unique(NodesEnriched);    % delete the Node in the Standard Finite Element Method
end
