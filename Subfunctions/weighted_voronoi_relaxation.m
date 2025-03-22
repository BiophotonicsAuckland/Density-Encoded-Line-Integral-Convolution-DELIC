function [points] = weighted_voronoi_relaxation(Initpoints, numIterations, overRelaxationFactor, scalarField, plotShow, verbose)
%% weighted_voronoi_relaxation
% 
% This function computed a centroidal voronoi tesselation using Lloyd's algorithm
%
% Inputs:
%   Initpoints           - Initial positions of the seed points. Each row
%                          corresponds to a point in the 2D domain.
%   numIterations        - (integer) Number of iterations of Lloyd's algorithm. 
%                          Default is 10 if not provided.
%   overRelaxationFactor - Does not work, just set to zero. Defaults to zero.
%   scalarField          - A 2D scalar field used to compute the weights for the 
%                          Voronoi cells. 
%   plotShow             - (binary) A flag to control whether to show intermediate plots during 
%                          the relaxation process (1 to plot, 0 to suppress).
%   verbose              - (binary) A flag to control whether to display iteration information 
%                          (1 to display, 0 to suppress).
%
% Outputs:
%   points               - (P x 2 matrix) The updated positions of the seed points after relaxation.
%                         
% Example:
%   % Example usage:
%   scalarField = rand(100, 100); % Random scalar field
%   Initpoints = [rand(100, 1) * 100, rand(100, 1) * 100]; % Random initial points
%   points = weighted_voronoi_relaxation(Initpoints, 100, 0.2, scalarField, 1, 1);

%%
    if nargin < 3
        numIterations = 10;
        overRelaxationFactor = 0;
    end

    if isempty(scalarField)
        scalarField = ones(size(max(Initpoints, 1), size(max(Initpoints, 2))));
    end

    % Generate meshgrid from input scalar field
    [xGrid, yGrid] = meshgrid(linspace(1, size(scalarField, 2), size(scalarField, 2)), linspace(1, size(scalarField, 1), size(scalarField, 1)));
    xLimits = [1, size(scalarField, 2)];
    yLimits = [1, size(scalarField, 1)];

    points = Initpoints;
    numPoints = size(points, 1);

    interpObj = griddedInterpolant(xGrid', yGrid', scalarField');
    epsilon = 1e-6;

    % Begin Voronoi Iteration
    for iter = 1:numIterations
        if verbose == 1
            disp(['Iteration number ', num2str(iter), ' of ', num2str(numIterations)])
        end

        % Generate Voronoi tessellation
        [V, C] = voronoin(points,{'Qz'});
        % {QZ} adds a point above the paraboloid of lifted sites for a Delaunay triangulation.
        % It allows the Delaunay triangulation of cospherical sites. It
        % reduces precision errors for nearly cospherical sites.
        
        newPoints = zeros(size(points));

        for i = 1:numPoints
            if all(C{i} ~= 1)  % If cell is bounded
                vertices = V(C{i}, :);
 
                % Clip cell to dom limit
                % vertices = clipVoronoiCell(vertices, xLimits, yLimits);
                
                if isempty(vertices)
                    newPoints(i, :) = points(i, :); % Reuse the old point if clip fail
                else
                    % Compute weighted centroid
                    weights = interpObj(vertices(:, 1), vertices(:, 2));
                    centroid = sum(vertices .* weights, 1) / sum(weights);
                    
                    % Over-relaxation step
                    points(i, :) = points(i, :) + overRelaxationFactor * (centroid - points(i, :));

                    % Check for NaN values
                    if any(isnan(centroid))
                        newPoints(i, :) = points(i, :); 
                    else
                        newPoints(i, :) = centroid;
                    end
                end
            else
                % Handle unbound cells 
                newPoints(i, :) = points(i, :);
            end
        end

        % Clip points to dom limits
        newPoints = min(max(newPoints, [xLimits(1), yLimits(1)]), [xLimits(2), yLimits(2)]);

        % Update and remove NaN and duplicate points
        points = unique(newPoints, 'rows');
        points = points(~any(isnan(points), 2), :);

        % Regenerate clipped points randomly across domain
        if size(points, 1) < numPoints
            numNewPoints = numPoints - size(points, 1);
            
            %Regenerate seeds randomly and scale to domain size
            x = [];y = [];
            timeout = 0;
            while numel(x) < numNewPoints
                timeout = timeout + 1 ;
                extra_points = numNewPoints - numel(x);
                extra_x = (xLimits(2) - xLimits(1)) * rand(extra_points, 1) + xLimits(1);
                extra_y = (yLimits(2) - yLimits(1)) * rand(extra_points, 1) + yLimits(1);

                extra_x = unique(round(extra_x / epsilon) * epsilon);
                extra_y = unique(round(extra_y / epsilon) * epsilon);

                x = [x; extra_x];
                y = [y; extra_y];

                % Remove duplicates 
                x = unique(round(x / epsilon) * epsilon);
                y = unique(round(y / epsilon) * epsilon);
                if timeout > 1000 %Force Exit after 1000 iterations
                    disp('Voronoi Relaxation Point Regenration Failed')
                    break
                end
            end

            newPoints = [x,y];
            points = [points; newPoints(1:numNewPoints, :)];

        end

        if plotShow == 1
            figure(1);
            plot(points(:, 1), points(:, 2), 'k.'); 
            hold on;
            [vx, vy] = voronoi(points(:, 1), points(:, 2));
            plot(vx, vy, 'r-', 'LineWidth', 1); 
            xlim(xLimits);
            ylim(yLimits);
            axis equal;
            title(['Voronoi Relaxation step num: ', num2str(iter)]);
            hold off;
        end
    end

    % Function to clip Voronoi cells to domain limits
    function clippedVertices = clipVoronoiCell(vertices, xLimits, yLimits)
        k = convhull(vertices);
        polyin = polyshape(vertices(k, 1), vertices(k, 2));
        rect = polyshape([xLimits(1) xLimits(1) xLimits(2) xLimits(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)]);
        polyout = intersect(polyin, rect);
        clippedVertices = polyout.Vertices;
    end
end
