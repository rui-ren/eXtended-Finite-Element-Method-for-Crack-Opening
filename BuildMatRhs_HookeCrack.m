% Hooke Srain Matrix 
function [MAT, RHS] = BuildMatRhs_HookeCrack(...
     NMat,  dNdxMat,  dNdyMat,  MMat,  dMdxMat,  dMdyMat, ...
    F1Mat, dF1dxMat, dF1dyMat, F2Mat, dF2dxMat, dF2dyMat, ...
    F3Mat, dF3dxMat, dF3dyMat, F4Mat, dF4dxMat, dF4dyMat, ...
    xxInt, yyInt, wwInt, Nodes, lambda, mu, fx, fy, nQ, NodeNum, MAT, RHS)

% Integrate the weak form in the domain, corresponding to a Hooke solid.

% Initialization.
ElemMAT11 = zeros(24, 24);
ElemMAT12 = zeros(24, 24);
ElemMAT21 = zeros(24, 24);
ElemMAT22 = zeros(24, 24);
ElemRHS1  = zeros(24, 1);
ElemRHS2  = zeros(24, 1);

% Loop over integration points.
for i = 1 : nQ

    N  =  NMat(:, i); Nx  =  dNdxMat(:, i); Ny  =  dNdyMat(:, i);
    M  =  MMat(:, i); Mx  =  dMdxMat(:, i); My  =  dMdyMat(:, i);
    F1 = F1Mat(:, i); F1x = dF1dxMat(:, i); F1y = dF1dyMat(:, i);
    F2 = F2Mat(:, i); F2x = dF2dxMat(:, i); F2y = dF2dyMat(:, i);
    F3 = F3Mat(:, i); F3x = dF3dxMat(:, i); F3y = dF3dyMat(:, i);
    F4 = F4Mat(:, i); F4x = dF4dxMat(:, i); F4y = dF4dyMat(:, i);

    NxNxT = [Nx; Mx; F1x; F2x; F3x; F4x] * [Nx; Mx; F1x; F2x; F3x; F4x]';    % The overall shape function.
    NxNyT = [Nx; Mx; F1x; F2x; F3x; F4x] * [Ny; My; F1y; F2y; F3y; F4y]';    % The overall shape function.
    NyNxT = [Ny; My; F1y; F2y; F3y; F4y] * [Nx; Mx; F1x; F2x; F3x; F4x]';    % The overall shape function.
    NyNyT = [Ny; My; F1y; F2y; F3y; F4y] * [Ny; My; F1y; F2y; F3y; F4y]';    % The overall shape function.

    % Compute element matrices.
    ElemMAT11 = ElemMAT11 + wwInt(i) * ( (lambda+2*mu) * NxNxT + mu * NyNyT );
    ElemMAT12 = ElemMAT12 + wwInt(i) * ( lambda * NxNyT + mu  * NyNxT );
    ElemMAT21 = ElemMAT21 + wwInt(i) * ( lambda * NyNxT + mu  * NxNyT );
    ElemMAT22 = ElemMAT22 + wwInt(i) * ( (lambda+2*mu) * NyNyT + mu * NxNxT );

    % Compute right hand side.
    ElemRHS1 = ElemRHS1 + wwInt(i) * ( [N; M; F1; F2; F3; F4] * fx );
    ElemRHS2 = ElemRHS2 + wwInt(i) * ( [N; M; F1; F2; F3; F4] * fy );

end

% Add element contribution to global matrix.
uuNodes = [Nodes          Nodes+2*NodeNum  Nodes+4*NodeNum  Nodes+5*NodeNum  Nodes+ 6*NodeNum  Nodes+ 7*NodeNum];
vvNodes = [Nodes+NodeNum  Nodes+3*NodeNum  Nodes+8*NodeNum  Nodes+9*NodeNum  Nodes+10*NodeNum  Nodes+11*NodeNum];

MAT(uuNodes, uuNodes) = MAT(uuNodes, uuNodes) + ElemMAT11;
MAT(uuNodes, vvNodes) = MAT(uuNodes, vvNodes) + ElemMAT12;
MAT(vvNodes, uuNodes) = MAT(vvNodes, uuNodes) + ElemMAT21;
MAT(vvNodes, vvNodes) = MAT(vvNodes, vvNodes) + ElemMAT22;

RHS(uuNodes) = RHS(uuNodes) + ElemRHS1;
RHS(vvNodes) = RHS(vvNodes) + ElemRHS2;
