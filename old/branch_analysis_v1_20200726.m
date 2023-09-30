% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% Segment & analyze dendritic y-branches %%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20200726: nks @ kwanlab, ysm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %gcp;

% Set paths
% addpath(genpath('/Users/neilsavalia/Documents/yale/kwanlab/project/dendrites_ketamine/scripts'));
addpath(genpath('/Users/neilsavalia/Documents/matlab'));

% % setup figure properties
% setup_figprop;

%% File loader

% get file
img_in = [];
while isempty(img_in)
    fl = imgetfile;
    [ff, fn, fp] = fileparts(fl);
    img_in = ScanImageTiffReader.ScanImageTiffReader(fl);
end

% preprocess
img = img_in.data;  % img is flipped L->R than rotated 90 from raw, readjusted in next line.
img = rot90(img(:, size(img, 2):-1:1, :));
hdrs = cellstr(splitlines(img_in.descriptions));
offsets = str2num(cell2mat(strrep(hdrs(contains(hdrs, 'c0=')), 'c0=', '')));
for i = 1:size(img, 3)
    img(:, :, i) = img(:, :, i) + offsets(i);
end

%% Plot & select

fig = figure('Position', [0, 0, 2000, 1000]); 

% full img max projection
s1 = subplot(2, 3, 1);
img_mean = mean(img, 3);
img_max = max(img, [], 3);
img_proj = img_max;
imagesc(img_proj);
axis off
colormap([zeros(500, 1), linspace(0, 1, 500)', zeros(500, 1)]);
colormapeditor;
title({'Image 3D Max Projection'; 'Sub-select X-Y (draw box) around branch then press ENTER'});

subsel = [];
branches = [];
while isempty(subsel)

    % draw an roi
    img_ss = drawrectangle;         % box will show up blue
    img_ss.Color = [.9, .3, .1];    % set completed ROIs to red/orange color

    % edit; do not move on until enter key is pressed
    input(['Setting down box containing Y-shaped branch... Press ''ENTER'' to proceed...'], 's');

    % get current roi mask & centroid
    disp(['Extracting data in bounding box...']);
    branches.subselected.vertices_x = round(img_ss.Vertices(:, 1));
    branches.subselected.vertices_y = round(img_ss.Vertices(:, 2));
    branches.subselected.boxmask = poly2mask(branches.subselected.vertices_x, branches.subselected.vertices_y, size(img, 1), size(img, 2));
    branches.subselected.boxproj = img_proj(branches.subselected.vertices_y(1):branches.subselected.vertices_y(2), branches.subselected.vertices_x(1):branches.subselected.vertices_x(3));
    branches.subselected.boxdata = img(branches.subselected.vertices_y(1):branches.subselected.vertices_y(2), branches.subselected.vertices_x(1):branches.subselected.vertices_x(3), :);
    
    % exit condition for while
    subsel = 1;

end

% subselected x & y
s2 = subplot(2, 3, 2);
imagesc(branches.subselected.boxproj);
axis off
colormap([zeros(500, 1), linspace(0, 1, 500)', zeros(500, 1)]);
title('Sub-selected Branch 3D Max Projection');

% branch mask ----> FINAL FORM SHOULD MAKE THIS WITH A THRESHOLD SLIDER THAT DROPS OR KEEPS ADDITIONAL PIXELS
branches.subselected.boxthresh = 27.4/100;   % threshold for subselected X-Y branch image after smoothing (used to grab branch ROI from bwlabel output
branches.subselected.boxprojsmooth = imgaussfilt(branches.subselected.boxproj, 2, 'FilterSize', 3, 'Padding', 'replicate', 'FilterDomain', 'spatial');  % smooth sub-selected branch
branches.fullbranch.bw_extract = bwlabel(branches.subselected.boxprojsmooth > (branches.subselected.boxthresh * max(max(img_proj))));  % extract connected components greater than threshold 
branches.fullbranch.bw_extract_tab = tabulate(branches.fullbranch.bw_extract(:));     % tabulate size of connected components
branches.fullbranch.branchmask = (branches.fullbranch.bw_extract == (find(branches.fullbranch.bw_extract_tab(:, 2) == max(branches.fullbranch.bw_extract_tab(2:end, 2))) - 1));     % take largest connected component (assumed to be the branch)
s3 = subplot(2, 3, 3);
imagesc(branches.fullbranch.branchmask);
axis off
colormap([zeros(500, 1), linspace(0, 1, 500)', zeros(500, 1)]);
title({'Sub-selected Branch; Binarized Mask'; 'Draw polygon around branch point then press ENTER'});

% select branch point to separate shafts
bpoint = [];
while isempty(bpoint)

    % draw an roi
    img_ss = drawpolygon;            % circle will show up blue
    img_ss.Color = [.9, .3, .1];    % set completed ROIs to red/orange color

    % edit; do not move on until enter key is pressed
    input(['Setting down polygon containing branchpoint... Press ''ENTER'' to proceed...'], 's');

    % get current roi mask & centroid
    disp(['Extracting data in bounding box...']);
    branches.fullbranch.branchpointx = img_ss.Position(:, 1);
    branches.fullbranch.branchpointy = img_ss.Position(:, 2);
    branches.fullbranch.branchpoint_center = [mean(branches.fullbranch.branchpointx), mean(branches.fullbranch.branchpointy)]; 
    branches.fullbranch.branchpointmask = poly2mask(branches.fullbranch.branchpointx, branches.fullbranch.branchpointy, size(branches.fullbranch.branchmask, 1), size(branches.fullbranch.branchmask, 2));
    
    % separate branches
    branches.compare.masks_all = bwlabel(branches.fullbranch.branchmask .* ~branches.fullbranch.branchpointmask);
    
    % exit condition for while
    bpoint = 1;

end

% fit splines to each separated branch towards branchpoint
s4 = subplot(2, 3, 4);
centroid = branches.fullbranch.branchpoint_center;
for i = 1:size(unique(branches.compare.masks_all(branches.compare.masks_all > 0)), 1)
    
    % current mask name (temp until angles can be compared)
    bname = eval(['''branch' num2str(i) '''']);
    bmask = (branches.compare.masks_all == i);
    
    % get coordinates in box and fit split, append branch point centroid
    [yc, xc]    = find(branches.compare.masks_all == i);     % all points within branch
    xc          = [xc; centroid(1)];
    yc          = [yc; centroid(2)];
    
    % fit spline to x-y data for branch, x-natural
    pp = polyfix(xc, yc, 3, xc(end), yc(end));
    py = polyval(pp, xc);
    spline_fit = sortrows([xc, py]);
    
    % remove xc & yc values crossing threshold of branchpoint
    if median(xc) < centroid(1)
        if median(yc) < centroid(2)
            spline_fit = spline_fit(round(spline_fit(:, 1), 4) <= round(centroid(1), 4) & round(spline_fit(:, 2), 4) <= round(centroid(2), 4), :);
        else
            spline_fit = spline_fit(round(spline_fit(:, 1), 4) <= round(centroid(1), 4) & round(spline_fit(:, 2), 4) >= round(centroid(2), 4), :);
        end
    else
        if median(yc) < centroid(2)
            spline_fit = spline_fit(round(spline_fit(:, 1), 4) >= round(centroid(1), 4) & round(spline_fit(:, 2), 4) <= round(centroid(2), 4), :);
        else
            spline_fit = spline_fit(round(spline_fit(:, 1), 4) >= round(centroid(1), 4) & round(spline_fit(:, 2), 4) >= round(centroid(2), 4), :);
        end
    end
    
    % plot branches with splines
    mm = imagesc(bmask .* i); hold on;
    set(mm, 'AlphaData', logical(bmask).*(4/5));
    plot(spline_fit(:, 1), spline_fit(:, 2), 'k', 'LineWidth', 2);
    
    % sink values
    eval(['branches.compare.' bname '.mask = bmask;']);
    eval(['branches.compare.' bname '.spline_fit = spline_fit;']);
    
    % clean up
    clear bname bmask spline_fit yc xc pp py mm
end
scatter(branches.fullbranch.branchpoint_center(1), branches.fullbranch.branchpoint_center(2), 200, 'r', 'filled');

% calculate angle between each spline, determine proximal (relative to soma) vs. two distal components
fnames = fieldnames(branches.compare);
bnames = fnames(contains(fnames, 'branch'));
x0 = branches.fullbranch.branchpoint_center(:, 1);
y0 = branches.fullbranch.branchpoint_center(:, 2);
offset = 4; % points to offset from branch point -- could use zero but slightly odd behavior for splines that converge on branch point from far away
for i = 1:size(bnames, 1)
    
    % get last point of spline (separated by offset) leading up to branch point center
    nb1 = setdiff(eval(['branches.compare.' bnames{i} '.spline_fit']), branches.fullbranch.branchpoint_center, 'rows');
    if pdist([eval(['branches.compare.' bnames{i} '.spline_fit(end, :)'])', branches.fullbranch.branchpoint_center']') < pdist([eval(['branches.compare.' bnames{i} '.spline_fit(1, :)'])', branches.fullbranch.branchpoint_center']')
        x1 = nb1(end - offset, 1);
        y1 = nb1(end - offset, 2);
    else
        x1 = nb1(1 + offset, 1);
        y1 = nb1(1 + offset, 2);
    end
    
    for j = 1:size(bnames, 1)
        
        % get last point of spline (separated by offset) leading up to branch point center
        nb2 = setdiff(eval(['branches.compare.' bnames{j} '.spline_fit']), branches.fullbranch.branchpoint_center, 'rows');
        if pdist([eval(['branches.compare.' bnames{j} '.spline_fit(end, :)'])', branches.fullbranch.branchpoint_center']') < pdist([eval(['branches.compare.' bnames{j} '.spline_fit(1, :)'])', branches.fullbranch.branchpoint_center']')
            x2 = nb2(end - offset, 1);
            y2 = nb2(end - offset, 2);
        else
            x2 = nb2(1 + offset, 1);
            y2 = nb2(1 + offset, 2);
        end
        
        % calculate angle
        v1 = [x1, y1, 0] - [x0, y0, 0];
        v2 = [x2, y2, 0] - [x0, y0, 0];
        branches.compare.thetas(i, j) = atan2(norm(cross(v1, v2)), dot(v1, v2));
        branches.compare.angles(i, j) = rad2deg(branches.compare.thetas(i, j));
        
        % write angle onto plot between mask centroids
        if i ~= j
            midangpt = mean([x1, y1; x2, y2], 1);
            cent = branches.fullbranch.branchpoint_center;
            dd = pdist([midangpt; cent]);
            tx = cent(1) - ((20 * (cent(1) - midangpt(1))) / dd);
            ty = cent(2) - ((20 * (cent(2) - midangpt(2))) / dd);
            text(tx, ty, [sprintf('%.2f', branches.compare.angles(i, j)) '' char(176)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            clear midangpt cent
        end
        
    end
    
end
clear x0 x1 x2 y0 y1 y2 offset

% re-attribute names based on angles (proximal, distal1, distal2)

% ROI timecourses proximal & distal branches
figure; 
for i = 1:3
    subplot(3, 1, i);
    tcr = reshape(box_tc, [size(box_tc, 1) * size(box_tc, 2), size(box_tc, 3)]); 
    dbr = reshape((branches.compare.masks_all == i), [size(branches.compare.masks_all, 1) * size(branches.compare.masks_all, 2), 1]);
    plot(mean(tcr(dbr, :), 1), 'k');
    if i == 1
        legend('proximal', 'boxoff');
    elseif i == 2
        legend('distal1', 'boxoff');
    elseif i == 3
        legend('distal2', 'boxoff');
    end
    ylim([0, 2000]);
    ylabel('Raw F'); %xlabel('Frames');
end


% make hull ROI around branches w/ sequential components, plot




