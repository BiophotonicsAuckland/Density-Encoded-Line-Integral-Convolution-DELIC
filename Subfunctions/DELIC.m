function [licImage,rc,texture] = DELIC(U,V,scalarField,initInterval,voronoiIteration,numSteps,stepSize,sigma,LB,Parallel,verbose)
    %% DELICT: Density Encoded Line Integral Convolution for visualizing vector fields

    % This function performs Density Encoded Line Integral Convolution (DELICT) to visualize 
    % optical axis in birefringent biological tissue.
    % Citation: [...]

    % Input arguments:
    %   U               - 2D array representing the x-component of the vector field (optical axis).
    %   V               - 2D array representing the y-component of the vector field (optical axis).
    %   scalarField     - 2D array representing the scalar field (e.g., birefringence or another property).
    %   initInterval    - Integer defining the initial interval for setting up the texture (default: 4).
    %   voronoiIteration- Integer defining the number of iterations for Voronoi relaxation (default: 20).
    %   numSteps        - Integer defining the number of steps for each streamline (default: 50).
    %   stepSize        - Scalar defining the step size for the integration (default: 0.1).
    %   sigma           - Scalar defining the standard deviation for the Gaussian window (default: 1).
    %   LB              - Scalar defining the lower bound for Sobol sequence scaling (default: 0.35).
    %   Parallel        - Integer flag (0 or 1) to toggle use of parfor or for loop in LIC compute (default: 0).
    %   verbose         - Integer flag (0 or 1) to control verbosity of the function (default: 1).
    %
    % Output:
    %   licImage        - The final output DELIC image [without any post-processing]
    %   rc              - Seed points after N iterations of Lloyd's algorithm
    %   texture         - Texture image which was colvolved with the optical axis field

    % Example usage:
    %   [licImage] = DELIC(U, V, scalarField, 4, 20, 50, 0.1, 1, 0.35, 1);
    
    %% Check inputs and for sane values
    if nargin < 3
        error('Not enough inputs. Atleast U and V component of optical axis and encoding field needed')
    end

    if nargin < 4, initInterval = 4; end
    if nargin < 5, voronoiIteration = 20; end
    if nargin < 6, numSteps = 50; end
    if nargin < 7, stepSize = 0.1; end
    if nargin < 8, sigma = 1; end
    if nargin < 9, LB = 0.35; end
    if nargin < 10, Parallel = 0; end
    if nargin < 11, verbose = 1; end

    assert(isnumeric(U) && ismatrix(U), 'U should be a numeric matrix');
    assert(isnumeric(V) && ismatrix(V), 'V should be a numeric matrix');
    assert(isnumeric(scalarField) && ismatrix(scalarField), 'scalarField should be a numeric matrix');
    assert(isequal(size(U), size(V)), 'U and V should have the same dimensions');
    assert(isequal(size(U), size(scalarField)), 'U, V and scalarField should have the same dimensions');
    assert(initInterval > 0, 'initInterval should be a positive integer');
    assert(voronoiIteration > 0, 'voronoiIteration should be a positive integer');
    assert(numSteps > 0, 'numSteps should be a positive integer');
    assert(stepSize > 0, 'stepSize should be a positive number');
    assert(sigma > 0, 'sigma should be a positive number');
    assert(LB >= 0 && LB <= 1, 'LB should be between 0 and 1');

    U = double(U);
    V = double(V);
    scalarField = double(scalarField);

    %% Identify  domain
    NX = size(U,2); % Size of the domain [Num Columns]
    NY = size(U,1); % Size of the domain [Num Rows]

    %%  Normalise scalar field if not yet done
    vmin = min(scalarField(:)); vmax = max(scalarField(:));
    scalarField = (scalarField - vmin) / (vmax -vmin);
    
    %% Generate initial seeding field
    % Parameters
    interval = initInterval; % Interval for setting texture values
    texture = zeros(NY, NX);
    
    for i = 1:interval:NX %iterate through columns 
        for j = 1:interval:NY %iterate through rows
           
            newX = i;
            newY = j;
    
            newX = max(min(newX, NX), 1);
            newY = max(min(newY, NY), 1);
    
            texture(newY, newX) = 1; %rand();
        end
    end
    
    textureBin = rand(NY, NX) < scalarField;
    textureBin = texture.*textureBin;
    textureBin(scalarField==0) = 0; %Remove any left over non-birefringent seeds
    
    %% Padd initial seeding field and refine using Lloyd's algorithm
    
    % Make biased resampled meshgrid
    [rows, cols] = find(textureBin);
    InitPoints = [cols, rows];

    %Pad seeding grid to set domain limits to avoid unbounded points 
    [Padded_Grid,EdgePoints] = Grid_Padding(InitPoints,size(scalarField, 2),size(scalarField, 1),NY,NX);
    plotShow = 0;
    iteration = voronoiIteration; overRelax = 0;

    if verbose == 1
        tic
    end
    [rc] = weighted_voronoi_relaxation(Padded_Grid, iteration, overRelax,scalarField,plotShow,verbose);
    if verbose == 1
        disp(['Voronoi Relaxation computed in ',num2str(toc), ' seconds'])
    end

    % Remove padding
    rc = setdiff(rc, EdgePoints, 'rows');
   
%% Generate texture image using seeding field and assign random gray values using Sobol quaisirandom sequencing

    binaryImage = zeros(NY, NX);

    % Start Sobol sequence gray value generator
    numPoints = length(rc);

    p = sobolset(2, 'Skip', 1e3, 'Leap', 1e2); 
    p = scramble(p, 'MatousekAffineOwen');

    sobolPoints = net(p, numPoints);
    lowerBound = LB;
    sobolValues = lowerBound + (1 - lowerBound) * sobolPoints(:, 1); % Scale values to range [lowerBound, 1]

    for k = 1:numPoints
        % Get the x and y coordinates
        x = rc(k, 1);
        y = rc(k, 2);

        % Find the corresponding indices in the image
        i = round(x);
        j = round(y);

        % Set the corresponding pixel in the binary image to a Sobol sequence value
        binaryImage(j, i) = sobolValues(k);
    end


    %% Show final sparse texture. This texture is used for the LIC computation
  
    texture = binaryImage;
    if  verbose ==1
        figure(111);
        imshow(texture,[0 1]); colormap(gray)
    end

    %% Prepare for the LIC 

    % % Initialize the LIC image
    % licImage = zeros(NY, NX);
    
    % Normalize the vector field
    mag = sqrt(U.^2 + V.^2);
    U = U ./ (mag + eps); % Avoid division by zero
    V = V ./ (mag + eps);
    
    windo = gausswin(2*numSteps,sigma);
    windo = (windo(1:numSteps));
    
    %% Begin the LIC

    if verbose == 1
        tic
    end
    
    if Parallel == 1
        licImage = LIC_Parfor(U,V,texture,windo,numSteps,stepSize);
    else
        licImage = LIC(U,V,texture,windo,numSteps,stepSize);
    end

    if verbose == 1
        tic
    end
end