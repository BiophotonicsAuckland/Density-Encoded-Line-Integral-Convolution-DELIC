function [Padded_Grid,EdgePoints] = Grid_Padding(grid,DomEdgeX,DomEdgeY,NumPadX,NumPadY)
% This funcion eats a sparse grid and padds it at some user set frequency and
% user set edge to prepare grid edges for Voronoi tesselation

    xEdges = linspace(1, DomEdgeX, NumPadX);
    yEdges = linspace(1, DomEdgeY, NumPadY);
    
    bottomEdgeX = xEdges;
    bottomEdgeY = ones(size(xEdges));
    
    topEdgeX = xEdges;
    topEdgeY = (DomEdgeY) * ones(size(xEdges));
    
    leftEdgeX = ones(size(yEdges));
    leftEdgeY = yEdges;
    
    rightEdgeX = (DomEdgeX) * ones(size(yEdges));
    rightEdgeY = yEdges;
    
    edgePointsX = [bottomEdgeX, topEdgeX, leftEdgeX, rightEdgeX];
    edgePointsY = [bottomEdgeY, topEdgeY, leftEdgeY, rightEdgeY];
    
    EdgePoints = [edgePointsX(:), edgePointsY(:)];
    EdgePoints = unique(EdgePoints, 'rows');
    Padded_Grid = unique([grid; EdgePoints], 'rows');

end
