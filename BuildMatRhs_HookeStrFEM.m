%
%
% MATRIX Assemble for standard Finite Element Method

function [MAT, RHS] = BuildMatRhs_HookeStrFEM(NMat,dNdxMat,dNdyMat,...
    xxInt, yyInt, wwInt, Nodes, lambda, mu, fx, fy,nQ, NodeNum, MAT, RHS)

% Integrate the weak form in the domain, corresponding to a Hooke solid.

% Initialization.
ElemMAT11 = zeros(4,4);
ElemMAT12 = zeros(4,4);
ElemMAT21 = zeros(4,4);
ElemMAT22 = zeros(4,4);
ElemRHS1 = zeros(4,1);
ElemRHS2 = zeros(4,1);

% Loop over integration points

for i = 1 : nQ
    
    N = NMat(:, i);
    Nx  = dNdxMat(:, i);
    Ny = dNdyMat(:, i);
    
    NxNxT = Nx * Nx';
    NxNyT = Nx * Ny';
    NyNxT = Ny * Nx';
    NyNyT = Ny * Ny';
   
    % Compute element matrices.
    
    ElemMAT11 = ElemMAT11 + wwInt(i) * ((lambda + 2*mu)*NxNxT+mu*NyNyT);
    ElemMAT12 = ElemMAT12 + wwInt(i) * (lambda * NxNyT + mu* NyNxT);
    ElemMAT21 = ElemMAT21 + wwInt(i) * (lambda * NyNxT + mu * NxNyT);
    ElemMAT22 = ElemMAT22 + wwInt(i) * (lambda + 2*mu) * NyNyT + mu * NxNxT;
    
     % Compute right hand side.
    ElemRHS1 = ElemRHS1 + wwInt(i) * ( N * fx );
    ElemRHS2 = ElemRHS2 + wwInt(i) * ( N * fy );
end

    
    uuNodes = [Nodes   ];
    vvNodes = [Nodes + NodeNum];
    
    MAT(uuNodes, uuNodes) = MAT (uuNodes, uuNodes) + ElemMAT11;
    MAT(uuNodes, vvNodes) = MAT(uuNodes, vvNodes) + ElemMAT12;
    MAT(vvNodes, uuNodes) = MAT(vvNodes, uuNodes) + ElemMAT21;
    MAT(vvNodes, vvNodes) = MAT(vvNodes, vvNodes) + ElemMAT22;
    
    RHS(uuNodes) = RHS(uuNodes) + ElemRHS1;
    RHS(vvNodes) = RHS(vvNodes) + ElemRHS2;




end
