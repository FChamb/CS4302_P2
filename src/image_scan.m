% ============================
% CS4302 Practical 2 - Images
% File: image_scan.m
% ============================
% Steps:
% 1) Read each nuts image.
% 2) Segment nuts from white background using grayscale thresholding.
% 3) Clean mask, then split touching nuts with distance transform + watershed.
% 4) Count total nuts and classify each nut (almond / walnut / cashew).
% 5) Save overlay and a CSV summary.

clear; close all; clc;

imgNames = { '../images/Easy.jpeg', '../images/Medium.jpeg', '../images/Hard.jpeg', '../images/Very_hard.jpeg', '../images/Extreme.jpeg' };
outDir   = '../out';
if ~exist(outDir, 'dir'), mkdir(outDir); end

allRows = {};   % rows for final CSV table

for i = 1:numel(imgNames)
    inName = imgNames{i};

    if isempty(inName)
        fprintf('No file found for image set %d, skipping.\n', i);
        continue;
    end

    fprintf('Processing %s ...\n', inName);

    % ----- read image and optionally downscale -----
    I = imread(inName);
    maxDim = 1600;
    scale = min(1, max(size(I,1), size(I,2)) / maxDim);
    if scale > 1
        % if image is larger than maxDim, shrink it
        I = imresize(I, 1/scale);
    end

    %% ===== SEGMENTATION (Saturation only) =====
    % Nuts are coloured; background is low-saturation.

    Ihsv = rgb2hsv(I);
    S = Ihsv(:,:,2);

    % simple threshold on saturation
    mask = S > 0.20;

    % clean-up: fill holes, smooth edges, remove tiny specks
    mask = imfill(mask, 'holes');
    mask = imopen(mask, strel('disk', 3));    % smooth small bumps
    mask = imclose(mask, strel('disk', 2));   % close small gaps
    mask = bwareaopen(mask, 200);            % drop very small blobs

    % final segmentation mask (no watershed)
    seg = imopen(mask, strel('disk', 1));
    seg = bwareaopen(seg, 200);

    %% ===== LABEL OBJECTS AND MEASURE FEATURES =====
    CC = bwconncomp(seg);
    fprintf('  Found %d raw regions (before size filtering).\n', CC.NumObjects);

    if CC.NumObjects == 0
        fprintf('  WARNING: no objects found, skipping this image.\n');
        continue;
    end

    % --- area-based filtering to remove small fragments ---
    areas   = cellfun(@numel, CC.PixelIdxList);
    medArea = median(areas);

    % keep only regions that are at least half the median area
    keepIdx = find(areas > 0.5 * medArea);

    keepMask = false(size(seg));
    for k = 1:numel(keepIdx)
        keepMask(CC.PixelIdxList{keepIdx(k)}) = true;
    end
    seg = keepMask;

    % recompute connected components AFTER filtering
    CC = bwconncomp(seg);
    fprintf('  After size filtering: %d candidate nuts.\n', CC.NumObjects);

    if CC.NumObjects == 0
        fprintf('  WARNING: no objects left after filtering, skipping.\n');
        continue;
    end

    props = regionprops(CC, ...
        'Area','Perimeter','Eccentricity','Solidity', ...
        'MajorAxisLength','MinorAxisLength','BoundingBox','Centroid');

    % texture cue (walnuts are rougher)
    Igray = rgb2gray(I);
    E = entropyfilt(Igray, true(9));           % local entropy
    statsEntropy = regionprops(CC, E, 'MeanIntensity');

    %% ===== CLASSIFICATION =====
    areas2      = [props.Area];
    medAreaObj  = median(areas2);              % approx almond size

    labels = strings(CC.NumObjects,1);

    for n = 1:CC.NumObjects
        Areg  = props(n).Area;
        P     = props(n).Perimeter + eps;
        ecc   = props(n).Eccentricity;
        sol   = props(n).Solidity;
        maj   = props(n).MajorAxisLength;
        minr  = props(n).MinorAxisLength + eps;
        ar    = maj / minr;                       % aspect ratio
        circ  = 4*pi*Areg/(P^2);                  % circularity
        ent   = statsEntropy(n).MeanIntensity;    % roughness
        normA = Areg / medAreaObj;                % area relative to almond size

        % ---- 1) Walnuts: clearly larger and rougher ----
        if (normA > 1.8 && ent > 4.5) || (normA > 2.2)
            labels(n) = "walnut";

        % ---- 2) Cashews: medium size, curved, lower circularity ----
        elseif (normA > 0.7 && normA < 1.8 && ...
                circ < 0.70 && sol < 0.96 && ecc > 0.60)
            labels(n) = "cashew";

        % ---- 3) Almonds: smaller, smooth, elongated, fairly oval ----
        elseif (ecc > 0.75 && sol > 0.90 && ar > 1.4 && circ >= 0.70)
            labels(n) = "almond";

        % ---- 4) Fallbacks: decide mainly by size & circularity ----
        else
            if normA >= 1.8 && ent > 4.3
                labels(n) = "walnut";
            elseif circ < 0.70
                labels(n) = "cashew";
            else
                labels(n) = "almond";
            end
        end
    end

    %% ===== VISUAL OUTPUTS =====
    % distance map (for the steps figure only)
    D = bwdist(~seg);

    % segmentation overlay
    overlay = labeloverlay(I, seg, 'Transparency', 0.6);
    base = erase(inName, {'.jpeg','.JPEG','.jpg','.JPG'});
    imwrite(overlay, fullfile(outDir, sprintf('%s_overlay.png', base)));

    % bounding boxes + class labels
    RGB = I;
    for n = 1:CC.NumObjects
        bb = props(n).BoundingBox;
        RGB = insertShape(RGB,'Rectangle',bb,'LineWidth',3);
        RGB = insertText(RGB, props(n).Centroid, labels(n), 'BoxOpacity',0.6);
    end
    imwrite(RGB, fullfile(outDir, sprintf('%s_labeled.png', base)));

    %% ===== CONSOLE SUMMARY =====
    totalNuts = CC.NumObjects;
    nAlmond   = sum(labels == "almond");
    nWalnut   = sum(labels == "walnut");
    nCashew   = sum(labels == "cashew");

    fprintf('  Total: %d | almonds: %d, walnuts: %d, cashews: %d\n\n', ...
        totalNuts, nAlmond, nWalnut, nCashew);

    allRows(end+1,:) = { inName, totalNuts, nAlmond, nWalnut, nCashew }; %#ok<SAGROW>

    %% ===== STEP-BY-STEP MONTAGE (for report) =====
    fig = figure('Name', inName, 'Color','w');
    tlo = tiledlayout(2,3, 'Padding','compact','TileSpacing','compact');

    nexttile; imshow(I);        title('Input');
    nexttile; imshow(S,[]);     title('HSV S channel');
    nexttile; imshow(mask);     title('Initial mask');
    nexttile; imshow(seg);      title('Final segments');
    nexttile; imshow(D,[]);     title('Distance map (no watershed)');
    nexttile; imshow(RGB);      title('Labeled nuts');

    exportgraphics(fig, fullfile(outDir, sprintf('%s_steps.png', base)));
    close(fig);
end

% ===== SAVE CSV SUMMARY =====
T = cell2table(allRows, ...
    'VariableNames', {'image','total','almond','walnut','cashew'});
writetable(T, fullfile(outDir,'nut_counts_summary.csv'));

disp('Done. Results saved in out/');
