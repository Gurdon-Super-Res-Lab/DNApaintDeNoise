% Example data analysis for Nature Communications manuscript titled
% "Utilizing the predictable binding kinetics of DNA-PAINT to denoise
% super-resolution images" by Sirinakis et al. 2025
clear variables; clc

%% load example data and assign variables

% load example data
load("ExampleData.mat")

% assign data columns to variables
X = Data(:, 1); % X position (pixels)
Y = Data(:, 2); % Y position (pixels)
T = Data(:, 3); % frame number
E = Data(:, 4); % localization precision (pixels)

%% load ThunderSTORM csv file

% % Uncomment this section to import data from a ThunderSTORM csv file.
% % Import is based on the csv format generated with the example dataset
% % published with the ThunderSTORM paper as of June 2025
% 
% % user select file to import
% [rFile, rPath] = uigetfile('*.csv');
% 
% % import data from ThunderSTORM csv file
% TSdata = readmatrix([rPath, rFile]);
% 
% % assign imported data to named variables
% X = TSdata(:, 2); % X position (pixels)
% Y = TSdata(:, 3); % Y position (pixels)
% T = TSdata(:, 1); % frame number
% E = TSdata(:, end)./100; % localization precision (pixels) assuming a 100 nm pixel size

%% 1. run dbscan clustering on example data

% set the search radius to use with dbscan
epsilon = median(E);

% set the number of points for dbscan clustering
minpts=5;

% run dbscan on the example data - expect this step to take a few minutes
% depending on your computer speed (it will be very slow - hours - if you
% use datasets with hundreds of thousands of x, y points)
fprintf(1, 'Running DBSCAN clustering...');
dbIdx = dbscan([X, Y], epsilon, minpts);
fprintf(1, 'DONE\n');

%% 2. second clustering step

% run the step2cluster function which takes each cluster found via dbscan
% (above) and looks for sub-clusters using gaussian template matching to
% find the number of clusters and their initial location (within the larger
% cluster) and passes that to the k-means clustering function. This
% partitions each db-cluster into single molecule sized sub-clusters.
fprintf(1, 'Running second level clustering...');
[gIdx, silh] = step2Cluster(X, Y, E, dbIdx);
fprintf(1, 'DONE\n');

%% link the blinks and do the Anderson-Darling Test

% pre-allocate an array to hold the p-value found by the anderson-darling
% test for each (x, y) point that is a member of a cluster found above
pVal = zeros(size(gIdx));

% set the maximum gap size (frames) for linking blinks in time
gapFrames = 5;

% set a threshold for the minimum number of binding events each cluster
% must contain
bindEventThres = 8;

% get the unique ID for each cluster found in the second clustering step
clustIDs = unique(gIdx);

% get the total number of clusters that will be tested
numClust = numel(clustIDs);

% loop over each cluster
fprintf(1, 'Evaluating each cluster...');
for i = 1:numClust

    % ensure we don't use the zero cluster ID which corresponds to
    % unclustered points
    if clustIDs(i) ~= 0

        % 
        linkIdx = linkBlinks(T(gIdx == clustIDs(i)), gapFrames);

        [dToffs, ~, ~, ~] = makeTs(T(gIdx == clustIDs(i)), linkIdx);

        % check that the current cluster has at least the minimum number of
        % binding events
        if numel(dToffs) >= bindEventThres - 1

            % run the anderson-darling test on the curret cluster
            % dark-times (dToffs) to test the null hypothesis that the data
            % Comes from the specified distribution, in this case
            % exponetial
            [~, clusterPval] = adtest(dToffs, 'Distribution', 'exp');
            
            % store the p-value for this cluster in the p-value results
            % array
            pVal(gIdx == clustIDs(i)) = clusterPval;

        end
    end
end
fprintf(1, 'DONE\n');

%% threshold and make figures

% set a threshold for the p-value. Clusters at or above this threshold are
% considered valid and accepted and clusters below this threshold are
% considered invalid and rejected
pValThres = 0.01;

% make a logical index for each (x, y) point that is true (1) if that point
% is part of an accepted cluster or false (0) if that point is not
pIDX = pVal >= pValThres;

% open a new figure window
figure

% set the figure window to be a tiled layout (multiple plots)
t = tiledlayout(1, 2);

% set the focus on the first tile and plot all the (x, y) positions
nexttile(t)
gscatter(X(:), Y(:), gIdx(:))
axis image
xticks([])
yticks([])
legend('off')
title('Raw Data')

% set the focus on the second tile and plot everything above the p-value
% threshold that is considered valid and accepted
nexttile(t)
gscatter(X(pIDX), Y(pIDX), gIdx(pIDX))
axis image
xticks([])
yticks([])
legend('off')
title('Denoised Data')

%% save data

save('Example Output.mat', 'pVal', 'gIdx', 'dbIdx', 'epsilon', 'gapFrames', 'pValThres', 'minpts')