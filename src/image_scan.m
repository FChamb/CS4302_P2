function result = segment_and_count(imgPath, outdir, tag)
    I = imread(imgPath);
    Iu = im2uint8(I);

    Ilab = rgb2lab(Iu);
    L = Ilab(:,:,1); a = Ilab(:,:,2); b = Ilab(:,:,3);

    chroma = hypot(a, b);
    score = rescale((100 - L)) .* rescale(chroma);

    level = graythresh(score);
    BW0 = imbinarize(score, level*0.9);

    BW = imopen(BW0, strel('disk', 3));
    BW = imclose(BW, strel('disk', 5));
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, 200);

    D = -bwdist(~BW);
    D(~BW) = -Inf;
    Lws = watershed(D);
    BWsep = BW; BWsep(Lws==0) = 0;

    CC = bwconncomp(BWsep);
    stats = regionprops(CC, 'Area','Eccentricity','Solidity','Perimeter','BoundingBox');
    count = CC.NumObjects;

    steps = {
        I,            'step1_rgb';
        mat2gray(score),'step2_score';
        BW0,          'step3_thresh';
        BW,           'step4_morph';
        label2rgb(Lws,'jet','k','shuffle'), 'step5_watershed';
        insertObjectAnnotation(I,'rectangle',vertcat(stats.BoundingBox), ...
            arrayfun(@(i) sprintf('%d',i), 1:numel(stats),'UniformOutput',false)), 'step6_labels'
    };
    for i=1:size(steps,1)
        imwrite(steps{i,1}, fullfile(outdir, sprintf('%s_%s.png', tag, steps{i,2})));
    end

    result.count = count;
    result.stats = stats;
    result.mask  = BWsep;
end

difficulties = {'Easy','Medium','Hard','Very_hard','Extreme'};
imgdir = '../images';
T = table('Size',[numel(difficulties) 2],'VariableTypes',{'string','double'}, ...
          'VariableNames',{'difficulty','total_nuts'});

for i=1:numel(difficulties)
    imgPath = fullfile(imgdir, [difficulties{i} '.jpeg']);
    tag = sprintf('seg_%s', difficulties{i});
    R = segment_and_count(imgPath, '../out', tag);
    T.difficulty(i) = difficulties{i};
    T.total_nuts(i) = R.count;
end

writetable(T, fullfile('../out','counts_total.csv'));
disp(T)

function [labels, labelRGB, counts] = classify_nuts(I, mask, stats)
    L = bwlabel(mask);
    n = numel(stats);
    labels = strings(n,1);

    Igray = rgb2gray(I);

    for i = 1:n
        S = stats(i);
        objMask = (L == i);

        A  = S.Area;
        P  = S.Perimeter + eps;
        ecc = S.Eccentricity;
        sol = S.Solidity;
        circVar = (P.^2) / (4*pi*A);

        bb = round(S.BoundingBox);
        bb(1:2) = max(bb(1:2),1);
        bb(3) = max(bb(3),1); bb(4) = max(bb(4),1);
        crop = Igray(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(3)-1);
        crop = im2uint8(mat2gray(crop));
        glcm = graycomatrix(crop, 'Offset',[0 1; -1 1; -1 0; -1 -1], 'Symmetric',true);
        props = graycoprops(glcm, {'Contrast','Homogeneity'});
        texContrast = mean(props.Contrast);
        texHomog    = mean(props.Homogeneity);

        if (ecc > 0.88) && (sol > 0.90) && (circVar > 1.2 && circVar < 2.2)
            labels(i) = "almond";
        elseif (sol < 0.90 && circVar >= 2.0) || (texHomog < 0.35 && ecc > 0.75)
            labels(i) = "cashew";
        else
            if (ecc < 0.85 && (texContrast > 3.0 || sol < 0.95))
                labels(i) = "walnut";
            else
                labels(i) = "walnut";
            end
        end
    end

    cmap = containers.Map( {'almond','walnut','cashew'}, ...
                           {uint8([255 140 0]), uint8([139 69 19]), uint8([255 215 0])}); % colours
    [H,W,~] = size(I);
    labelRGB = I;
    L = bwlabel(mask);
    for i = 1:n
        c = cmap(labels(i));
        labelRGB(repmat(L==i,1,1,3)) = reshape(c,1,1,3);
    end

    counts = struct('almond', sum(labels=="almond"), ...
                    'walnut', sum(labels=="walnut"), ...
                    'cashew', sum(labels=="cashew"));
end

T2 = table('Size',[numel(difficulties) 4], 'VariableTypes',{'string','double','double','double'}, ...
           'VariableNames',{'difficulty','almond','walnut','cashew'});

for i=1:numel(difficulties)
    imgPath = fullfile(imgdir, [difficulties{i} '.jpg']);
    R = segment_and_count(imgPath, '../out', sprintf('seg_%s', difficulties{i}));

    I = imread(imgPath);
    [labels, labelRGB, counts] = classify_nuts(I, R.mask, R.stats);
    imwrite(labelRGB, fullfile('../out', sprintf('classes_%s.png', difficulties{i})));

    T2.difficulty(i) = difficulties{i};
    T2.almond(i) = counts.almond;
    T2.walnut(i) = counts.walnut;
    T2.cashew(i) = counts.cashew;
end
writetable(T2, fullfile('../out','counts_by_type.csv'));
disp(T2)