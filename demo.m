%% Load data
close all
clearvars

dn = load('BirefringenceEncodingField.mat');
dn = dn.scalarField; % Normalised Birefringence scalar field [0:1]
OpAx = load('OpAx.mat');
OpAx = OpAx.thet; % Optical axis field [0:180 degrees]

figure(101);imagesc(OpAx);colormap('hsv');colorbar() %Use Crameri's romao cyclic colour map for better visualisation
figure(102);imagesc(dn);colormap('turbo');colorbar() %Use Ametrine colormap from MATLAB File Exchange

%% Convert optical axis orientation to X and Y components
U = cosd(OpAx);
V = -sind(OpAx);
%% DELIC

%Define param
initInterval = 4;
voronoiIteration = 20;
numSteps = 50;
stepSize = 0.1;
sigma = 1;
LB = 0.35;
Parallel = 1;
verbose = 1;

[licImage,~,~] = DELIC(U,V,dn,initInterval,voronoiIteration,numSteps,stepSize,sigma,LB,Parallel,verbose);

%% Normalise DELIC image 
Vmax = max(licImage(:));
Vmin = min(licImage(:));
licImage = licImage - Vmin;
licImage = licImage / Vmax;
licImage = imadjust(licImage);
figure(112);imshow(licImage); colormap(gray)

%% Optional: Apply colour to DELIC image 
figure(113);
cmap = colormap('hsv'); % We recommend Crameri's romao cyclic colour map OR one of the cyclic maps in colorcet 
colorImage = GenerateColorDELIC(OpAx, licImage, cmap);
imshow(imlocalbrighten(colorImage,0.15));


