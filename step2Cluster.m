function [gIdx, silh] = step2Cluster(X, Y, eps, dbIdx)
% This function performs a second layer of clustering based on k-means and
% morhological features of the DBSCAN clusters

% INPUT
% X: list of x positions of blinks
% Y: list of y positions of blinks
% eps: list of localization precision of blinks
% dbIdx: Cluster index from dbScan analysis

% OUTPUT
% gIdx: global unique index for all blinks after the 2nd clustering
% silh: silhoutte score

% ensure that the array sizes match
if size(X, 1) == size(dbIdx, 1)

    % initialize output arrays
    gIdx=zeros(size(dbIdx));
    
    % initialize (x, y) values for the template image for xcorr
    w = 2;
    x= -w:1:w;
    y = x;
    [Xt, Yt] = meshgrid(x, y);
    
    % set the zoom for cluster image reconstruction
    theZoom = 10;
    
    % initalize silhoutte score
    silh = nan(max(dbIdx), 1);

    % loop over each DBSCAN cluster index
    for i = 1:max(dbIdx)   

        % select data that belongs only to a specific DBSCAN cluster
        sIdx = i == dbIdx;

        % set max/min values for image reconstruction of the DBSCAN
        % clusters add some additional (w + 1) padding for the xcorr function
        tmpXmin = min(X(sIdx)) - (w+1)/theZoom;
        tmpXmax = max(X(sIdx)) + (w+1)/theZoom;
        tmpYmin = min(Y(sIdx)) - (w+1)/theZoom;
        tmpYmax = max(Y(sIdx)) + (w+1)/theZoom;

        % set the image reconstruction blur amount
        imBlur = theZoom*median(eps(sIdx));
        
        % setup bin edges
        xEdges = tmpXmin:1/theZoom:tmpXmax;
        yEdges = tmpYmin:1/theZoom:tmpYmax;

        % make a 2D histogram of the (x, y) points to reconstruct the image
        % of this cluster
        [tmpIm] = histcounts2(X(sIdx), Y(sIdx), xEdges, yEdges);

        % add some gaussian blur
        tmpIm = imgaussfilt(tmpIm, imBlur);
        
        try            
            % make 2D gaussian template
            tBlur = imBlur;
            imT = (1/(2*pi*tBlur*tBlur))*exp(-(((Xt.^2)/(2.*(tBlur)^2))+((Yt.^2)/(2.*(tBlur)^2))));
            
            % normalized cross-correlation of tempate and cluster image
            tmpCorr = normxcorr2(imT, tmpIm);

            % crop correlated image to match original
            tmpCorr = tmpCorr(1 + w:1:end - w,1 + w:1:end - w);
            
            % apply threshold to the normalized cross-correlation image at
            % 0.1 level
            tmpCorr = max(tmpCorr, 0.1);
      
            %--------------------------------------------------------------
            % find peaks in the cross-correlation image and use their
            % locations as starting points for k-means clustering

            % find local maxima and generate a binary image
            BW = imregionalmax(tmpCorr, 8);

            % clear the borders
            BW([1,2], :) = false;
            BW([end-1, end], :) = false;
            BW(:, [1, 2]) = false;
            BW(:, [end-1, end]) = false;            

            % find the number of peaks
            nP = sum(BW(:)); 

            % find the position of the peaks
            [xP, yP] = find(BW);
            xP = tmpXmin + xP./theZoom;
            yP = tmpYmin + yP./theZoom;

            %--------------------------------------------------------------
            % run kmeans clustering using the number of clusters and
            % initial guess locations found above
            
            [tmpIdx, ~] = kmeans([X(sIdx), Y(sIdx)], nP, 'Start', [xP,yP], 'Distance', 'sqeuclidean', 'MaxIter', 1000);            

            % update global unique index
            gIdx(sIdx) = max(gIdx) + tmpIdx;

            % calculate silhouette score for each (x, y) point
            tempSil = silhouette([X(sIdx), Y(sIdx)], tmpIdx);

            % store the mean silhouette score
            silh(i, 1) = mean(tempSil);


        catch
            
            fprintf('Error cross-correlating the cluster with the template image.\n')
        end       
        
       
    end 
        
else
    fprintf('Datasets do not match.\n')
end