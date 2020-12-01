function [P,E,T,gradquotient] = flipped_mesh(Nvoxels)
%FLIPPED_MESH creates a flipped mesh compared to basic_mesh(1, Nvoxels)
%   The mesh is rotated 90 degrees and has the diagonal line creating
%   the triangulation going from the upper left corner to the lower right
%   corner.
    model = createpde(1);
    [X,Y] = meshgrid(-1:2/(Nvoxels-1):1);
    X = X(:);
    Y = Y(:);
    nodes = [X';Y'];
    K = delaunayTriangulation(nodes');
    elements = K.ConnectivityList';
    geometryFromMesh(model,nodes,elements);
    [P,E,T] = meshToPet(model.Mesh);
    gradquotient = 1;
end