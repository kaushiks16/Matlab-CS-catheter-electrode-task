% Task 1
% Read the DICOM file
X = dicomread('high-dose.dcm'); % FIrst step reading the DICOM image sequence file
% It gives an array of X:Y:1:frame so we need to extract a singular image

frameIdx = 10; % Sets the frame number Im gonna see 
frame = squeeze(X(:,:,1,frameIdx));      % Extracts first channel of frame  10 and produces 2D image with matrix X:Y
frame = double(frame);                   % convert to double data type - easier to work with I believe 

figure; %Opens a new figure window for plotting
imshow(frame, []);                       % auto-scales the display
axis image off;                          % makes it square 
colormap(gray(256)); %makes it grayscale
colorbar; % Adds the colourbar so we can see the display intensity of pixels


% Task 2 + Task 3

function [rows, cols, im3] = detectBlobPeaksROI(frame, N, roi, minStrengthPct, minDist)
% We declare a function called detectBlobpeaks that returns
% rows,columns(the coordinates of peaks) and im3 (the Laplacian response for visualisation and debugging)
% We input frame which is image, N which is number of peaks, roi which is
% region of interest this is important to remove outliers, miStrengthPct
% also important in removing and outliers and finally minDist which is min
% separation

% so first as I was getting errors I have put default arguments
    if nargin < 5 || isempty(minDist), minDist = 8; end %If minDist isnt provided it sets default of 8
    if nargin < 4 || isempty(minStrengthPct), minStrengthPct = 50; end % If minstrengthpct isnt provided it defaults to 50

    % 1) Image Smoothing and Laplacian response
    f = fspecial('gaussian', 9, 3); % Kernerl to smooth image and reduce noise so it can detect the max points
    newImg = imfilter(frame, f, 'replicate'); % Applies Gaussian filter to frame
    h = fspecial('laplacian', 0.2); %Kernel to measure curvature
    im3 = imfilter(newImg, h, 'replicate'); %Applies Laplacian to smoothed image to produce im3

    % 2) Finding the regional maximums (local peaks)
    peakMask = imregionalmax(im3); 

    % 3) Apply ROI mask to peakMask (so we only consider peaks inside ROI)
    xmin = roi(1); xmax = roi(2); ymin = roi(3); ymax = roi(4);%Unpacks the ROI vector into variables
    [rAll, cAll] = find(peakMask);% Finds row and column indices of all detected regional maxima in whole image
    insideROI = cAll >= xmin & cAll <= xmax & rAll >= ymin & rAll <= ymax;%We are only using peaks whose column (cAll) and row (rAll) fall inside the ROI rectangle
    rAll = rAll(insideROI);%Keeps only the row indices inside the ROI.
    cAll = cAll(insideROI);%Keeps only the column indices inside the ROI.

    % 4) Get strengths for those peaks
    linearIdx = sub2ind(size(im3), rAll, cAll);%We convert the voordinate subscripts of the ROI peaks into linear indices for indexing into im3
    strengths = im3(linearIdx);%Extracts the Laplacian response values at each candidate peak location

    % 5) Strength threshold (percentile) to remove weak noise peaks
    if ~isempty(strengths) && ~isempty(minStrengthPct)%Checks that there are candidate strengths and percentile thresholds provided
        thr = prctile(strengths, minStrengthPct);%Computes the threshold value equal to the minStrengthPct percentile of the candidate strengths.
        keep = strengths >= thr;%Selects peaks whose strength is at or above percentile threshold
        rAll = rAll(keep);%keeps only row coordinates of peaks that pass strength threshold
        cAll = cAll(keep);%keeps only the column indices of peaks that pass the strength threshold
        strengths = strengths(keep);%Keeps only the corresponding strength values
    end

    % 6) Further removal of anomalies
    % Build an image with peak strengths at peak locations
    peakImg = zeros(size(im3));
    if ~isempty(rAll)%Checks whether any peaks remain after thresholding.
        peakImg(linearIdx(keep)) = strengths; % This line places strengths into peakImg using the original linearIdx and keep
        linearIdxKept = sub2ind(size(im3), rAll, cAll);%Recomputes linear indices for the kept peaks (correct indices after filtering).
        peakImg = zeros(size(im3));%Resets peakImg to zeros to ensure only the kept peaks are placed.
        peakImg(linearIdxKept) = strengths;%Places the kept peak strengths into peakImg at their linear indices.
    end

    if ~isempty(rAll)
        se = strel('disk', round(minDist));%Creates a structuring element (disk) with radius minDist (rounded). This defines the neighborhood used for non‑maximum suppression.
        localMaxMask = peakImg == imdilate(peakImg, se) & peakImg > 0;%Performs non‑maximum suppression: dilate peakImg with the disk and 
        % keep only pixels that equal the dilated image (local maxima within the neighborhood) and are positive. This removes nearby weaker peaks, leaving only the strongest within each neighborhood.
        [rKept, cKept] = find(localMaxMask);%Finding row/column indices of the final local maxima after suppression.
        sKept = peakImg(localMaxMask);%Extracts the strengths of those final kept peaks.
    else
        rKept = []; cKept = []; sKept = [];%If no peaks remain output as empty arrays.
    end

    % 7) Sort by strength and pick top N
    if ~isempty(sKept) %Checks whether any peaks survived suppression.
        [sSorted, idxSort] = sort(sKept, 'descend');%Sorts the kept strengths in descending order and returns sort indices.
        idxSort = idxSort(1:min(N, numel(idxSort)));%Truncates the sorted index list to the top N peaks (or fewer if fewer peaks exist).
        rows = rKept(idxSort);%Selects the row coordinates of the top N peaks.
        cols = cKept(idxSort);%Selects the column coordinates of the top N peaks.
    else
        rows = [];
        cols = [];
    end

end
%End of function now to call upon it

% Parameters
numFrames = size(X,4);% Determines the number of frames in the DICOM sequence
Npeaks = 10; %Maximum number of peaks to request per frame.
xmin = 220; xmax = 420;    % X bounds (cols)(Found this through total image size and where the points are about)
ymin = 290; ymax = 500;    % Y bounds (rows)
roi = [xmin xmax ymin ymax];%Packs the ROI into a vector for passing to the function.
minStrengthPct = 80;      % keep peaks stronger than the 80th percentile 
minDist = 7;             % minimum distance between peaks (pixels)so all 10 arent showing the same electrode
pauseTime = 0.25;         % pause between frames for visualization when presenting

allCoords = cell(numFrames,1);%Preallocates a cell array to store detected coordinates for each frame.

for fIdx = 1:numFrames%Loop over every frame index.
    frame = squeeze(X(:,:,1,fIdx));%Extract the current frame's first channel and squeeze to 2‑D.
    frame = double(frame);%Convert to double cause its easier to work with

    [rows, cols, im3] = detectBlobPeaksROI(frame, Npeaks, roi, minStrengthPct, minDist);%Call function we made

    coords = [cols, rows];   % [X Y] pairs of coordinates assembled
    allCoords{fIdx} = coords;%Stores coordinates for this frame in cell array

    % Print results
    fprintf('Frame %d: %d points inside bbox\n', fIdx, size(coords,1));%Print a summary line showing how many points found inside ROI for the frame
    if isempty(coords)% If no points it prints no points otherwise it displays coordinate list
        disp('  (no points)');
    else
        disp(coords);
    end

    % Visualize: show original frame with filtered points
    figure(1); imshow(frame, []); hold on;%Shows the current frame in figure 1 and hold on to overlay markers.
    plot(coords(:,1), coords(:,2), 'r+', 'MarkerSize', 12, 'LineWidth', 2);%Plot detected points as red plus signs at their X,Y positions.
    rectangle('Position',[xmin, ymin, xmax-xmin, ymax-ymin], 'EdgeColor','g', 'LineWidth',1);%Draw the ROI rectangle in green for visual reference.
    title(sprintf('Frame %d — Filtered Electrode Coordinates', fIdx));%Add a title showing the frame number.
    hold off;
    drawnow;
    pause(pauseTime);%Release the hold, force the figure to update immediately, and pauses briefly so we can view the frame.
end

%Task4

% Parameters 
numFrames = numel(allCoords);%Recompute number of frames from allCoords
markerSize = 40;%Marker size for scatter3 plotting.
maxLinkDist = 25;   % max distance (pixels) to link points between frames for trajectories
showTrajectories = true;  % Flag to enable plotting of linked trajectories.

% Build arrays for scatter3
X = []; Y = []; Z = [];%Initialize arrays to accumulate all detected X, Y and frame indices for scatter plotting.
for f = 1:numFrames %Loop to collect all detections into the X, Y, Z arrays. Z stores the frame index for each point.
    pts = allCoords{f};   % Nx2 [x y]
    if isempty(pts), continue; end
    X = [X; pts(:,1)];
    Y = [Y; pts(:,2)];
    Z = [Z; repmat(f, size(pts,1), 1)];
end

% 3D scatter: X (cols), Y (rows), Z (frame index)
figure;
scatter3(X, Y, Z, markerSize, Z, 'filled');%Create a 3D scatter plot where color is mapped to frame index Z. Points are filled circles.
colormap(jet);%Use the jet colormap and show a colorbar to indicate frame index mapping
colorbar;
xlabel('X (cols)');
ylabel('Y (rows)');
zlabel('Frame index (time)');
title('Electrode detections across frames (Z = frame index)');
axis tight;%Tighten axis limits to the data range.
set(gca,'ZDir','reverse');   %reverse Z so earlier frames appear "on top"
grid on;
view(45,30);%Turn on grid and set a 3D view angle.

% bounding box projection for reference (first frame)
hold on;
xmin = 250; xmax = 440; ymin = 300; ymax = 470;
% Prepare to overlay reference rectangles at z=0 and z=numFrames; set a slightly different ROI for display.
plot3([xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], zeros(1,5), 'g--', 'LineWidth', 1);
%Plot the ROI rectangle at z=0 as a dashed green line.
plot3([xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], numFrames*ones(1,5), 'g--', 'LineWidth', 1);
%Plot the same rectangle at z=numFrames for context.

% Trajectory linking (using greedy nearest-neighbor)
if showTrajectories
    % Initialize tracks with points from first non-empty frame
    tracks = {};   % each track: Nx3 [x y z]
    firstIdx = find(~cellfun(@isempty, allCoords), 1, 'first');%Initialize an empty cell array tracks. Find the first frame that has any detections.
    if ~isempty(firstIdx)%Only proceed if there is at least one non‑empty frame.
        pts = allCoords{firstIdx};%Initialize one track per detected point in the first non‑empty frame. Each track stores rows of [x y frame].
        for p = 1:size(pts,1)
            tracks{end+1} = [pts(p,1), pts(p,2), firstIdx]; 
        end

        % For each subsequent frame, link points to existing tracks
        for f = firstIdx+1:numFrames%Loop over subsequent frames to assign current points to existing tracks.
            pts = allCoords{f};
            if isempty(pts), continue; end%Skip frames with no detections.

            % Track last positions
            lastPos = cellfun(@(t) t(end,1:2), tracks, 'UniformOutput', false);%Extract the last known (x,y) position of every existing track and stack into a matrix for distance computation.
            lastPosMat = vertcat(lastPos{:});   % Tx

            % Compute distances between lastPos and current pts
            D = pdist2(lastPosMat, pts);   % Compute pairwise Euclidean distances between each track's last position and each current point (TxM matrix).

            % Greedy linking: for each track, find nearest point within maxLinkDist
            assignedPts = false(size(pts,1),1);%Logical array to mark which current points have been assigned to a track.
            %For each existing track, find the nearest current point. If the nearest point is within maxLinkDist and not already assigned, append it to the track; otherwise append a NaN row to keep time alignment (track had no assignment this frame).
            for t = 1:size(D,1)
                [dmin, idxMin] = min(D(t,:));
                if dmin <= maxLinkDist && ~assignedPts(idxMin)
                    % append to track
                    tracks{t}(end+1, :) = [pts(idxMin,1), pts(idxMin,2), f];
                    assignedPts(idxMin) = true;
                else
                    % no assignment: optionally append NaN row to keep time alignment
                    tracks{t}(end+1, :) = [NaN, NaN, f];
                end
            end

            % Any unassigned points -> start new tracks
            for p = 1:size(pts,1)
                if ~assignedPts(p)
                    tracks{end+1} = [pts(p,1), pts(p,2), f];
                end
            end
        end

        % Plot tracks as lines in 3D
        for t = 1:numel(tracks)
            tr = tracks{t};
            % After linking, plot each track as a 3D line, skipping tracks with fewer than two valid points. NaN rows break continuity.
            valid = ~isnan(tr(:,1));
            if sum(valid) >= 2
                plot3(tr(valid,1), tr(valid,2), tr(valid,3), '-','LineWidth',1.5);
            end
        end
    end
end

hold off;
%Release the plot hold


% Task 5

% Parameters (set according to the file you used)
pixel_size_mm = 0.25;    % mm per pixel (Table 1)
frame_rate = 7.5;        % frames per second for high-dose.dcm (change if using low-dose)
maxLinkDist = 25;        % pixels, max linking distance between frames (tune if needed)

% Build tracks by greedy nearest-neighbour linking
numFrames = numel(allCoords);%Recompute numFrames and reinitialize tracks for a second pass that builds tracks for statistical analysis.
tracks = {};   % each track is an Nx3 array [x y frame]

%1. Find the first non‑empty frame and error out if there are no detections at all.
firstIdx = find(~cellfun(@isempty, allCoords), 1, 'first');
if isempty(firstIdx)
    error('allCoords is empty: no detections found.');
end

% initialize tracks with first non-empty frame
pts = allCoords{firstIdx};   % Nx2 [x y]
for p = 1:size(pts,1)
    tracks{end+1} = [pts(p,1), pts(p,2), firstIdx]; %#ok<SAGROW>
end

% link subsequent frames
for f = firstIdx+1:numFrames
    pts = allCoords{f};   % Mx2
    if isempty(pts)
        for t = 1:numel(tracks)
            tracks{t}(end+1, :) = [NaN, NaN, f];
        end
        continue;
    end

    lastPos = cellfun(@(t) t(end,1:2), tracks, 'UniformOutput', false);
    lastPosMat = vertcat(lastPos{:});   % Tx2
    D = pdist2(lastPosMat, pts);        % TxM

    assignedPts = false(size(pts,1),1);
    for t = 1:size(D,1)
        [dmin, idxMin] = min(D(t,:));
        if dmin <= maxLinkDist && ~assignedPts(idxMin)
            tracks{t}(end+1, :) = [pts(idxMin,1), pts(idxMin,2), f];
            assignedPts(idxMin) = true;
        else
            tracks{t}(end+1, :) = [NaN, NaN, f];
        end
    end

    for p = 1:size(pts,1)
        if ~assignedPts(p)
            tracks{end+1} = [pts(p,1), pts(p,2), f];
        end
    end
end

% For each track compute how many valid (non‑NaN) samples it has and, if at least two, compute the mean X coordinate. Tracks with fewer than two valid samples keep -Inf for meanCols.
numTracks = numel(tracks);
meanCols = -Inf(numTracks,1);
validCounts = zeros(numTracks,1);
for t = 1:numTracks
    valid = ~isnan(tracks{t}(:,1)) & ~isnan(tracks{t}(:,2));
    validCounts(t) = sum(valid);
    if validCounts(t) >= 2
        meanCols(t) = mean(tracks{t}(valid,1));  % mean X
    end
end

%If no track has >=2 samples, fall back to track with most valid samples
if all(isinf(meanCols))
    [~, idxMax] = max(validCounts);
    tipIdx = idxMax;
else
    [~, tipIdx] = max(meanCols);
end
% diagnostics
fprintf('Chosen tip track index: %d (valid samples = %d)\n', tipIdx, validCounts(tipIdx));

% Extract the chosen tip track, select only valid rows, and separate X, Y and frame indices for statistical calculations.
tipTrack = tracks{tipIdx};        % Nx3 [x y frame]
valid = ~isnan(tipTrack(:,1)) & ~isnan(tipTrack(:,2));
tipX = tipTrack(valid,1);   % columns (X)
tipY = tipTrack(valid,2);   % rows    (Y)
tipFrames = tipTrack(valid,3);

% Compute mean and standard deviation of tip column (X) and row (Y) in pixels.
mean_col_px = mean(tipX);
std_col_px  = std(tipX);
mean_row_px = mean(tipY);
std_row_px  = std(tipY);
%Compute horizontal and vertical movement ranges in pixels.
horiz_range_px = max(tipX) - min(tipX);
vert_range_px  = max(tipY) - min(tipY);

% Convert means and standard deviations from pixels to millimetres using the pixel size.
mean_col_mm = mean_col_px * pixel_size_mm;
std_col_mm  = std_col_px  * pixel_size_mm;
mean_row_mm = mean_row_px * pixel_size_mm;
std_row_mm  = std_row_px  * pixel_size_mm;
%Convert movement ranges to millimetres.
horiz_range_mm = horiz_range_px * pixel_size_mm;
vert_range_mm  = vert_range_px * pixel_size_mm;

% --- Display results ---
fprintf('Tip electrode statistics (N = %d valid samples):\n', numel(tipX));
fprintf('  Mean column (X): %.2f px  (%.3f mm)\n', mean_col_px, mean_col_mm);
fprintf('  Std  column (X): %.2f px  (%.3f mm)\n', std_col_px,  std_col_mm);
fprintf('  Mean row    (Y): %.2f px  (%.3f mm)\n', mean_row_px, mean_row_mm);
fprintf('  Std  row    (Y): %.2f px  (%.3f mm)\n', std_row_px,  std_row_mm);
fprintf('\nMovement ranges:\n');
fprintf('  Horizontal range: %.2f px  (%.3f mm)\n', horiz_range_px, horiz_range_mm);
fprintf('  Vertical   range: %.2f px  (%.3f mm)\n', vert_range_px, vert_range_mm);

% Compute frame-to-frame speeds (mm/s) 
if numel(tipX) >= 2
    dx_px = diff(tipX);
    dy_px = diff(tipY);
    dist_px = sqrt(dx_px.^2 + dy_px.^2);
    dist_mm = dist_px * pixel_size_mm;
    dt = diff(tipFrames) / frame_rate;   % seconds between samples (handles skipped frames)
    velocity_mm_s = dist_mm ./ dt;
    fprintf('\nMean frame-to-frame speed: %.3f mm/s (median %.3f mm/s)\n', mean(velocity_mm_s,'omitnan'), median(velocity_mm_s,'omitnan'));
end
