function imagesc_ex(I, C, labx, laby, direction)

% plot imagesc image with x- and y-ticks on the right places. Provide the
% image to plot and the xand y labels for each of the axes (have to be the
% right length, matching the size of the image!)
% Note: xlabels and ylables are numerical. Ticks will be placed at logical
% values
% I:          2D double matrix
% C:          [low high] for plotting color range. Pass [] for automatic sclaing
% xlab, ylab: vectors of double with values of axes
% direction:  pass normal for normal plotting direction (image standard has
%             reversed y-dir)

if nargin<4
    error('Must pass image I, xlabels xlab and ylabels ylab. C can be passed as [] to plot without scale')
end

if length(labx) ~= size(I,2) || length(laby) ~= size(I,1)
    error('size of labels does not match up with the image')
end

if nargin<5
    direction = 'normal';
end

if isempty(C)
    imagesc(I);
else
    imagesc(I, C);
end

% calculate where ticks should be placed. Use below heuristic:
resx = 10^floor(log10(range(labx)));
resy = 10^floor(log10(range(laby)));
% determine the tick locations
atx = floor(min(labx)/resx)*resx:resx:10*resx;
aty = floor(min(laby)/resy)*resx:resy:10*resy;
% remove ticks outside range
atx = atx(atx>=min(labx) & atx<=max(labx));
aty = aty(aty>=min(laby) & aty<=max(laby));

% approximate to image label values
[~, ytick_idx] = min(abs(laby' - aty), [], 1);  % Find closest indices
[~, xtick_idx] = min(abs(labx' - atx), [], 1);

% Apply ticks and labels
set(gca, 'YTick', ytick_idx, 'YTickLabel', aty);
set(gca, 'XTick', xtick_idx, 'XTickLabel', atx);

% reverse direction of image to normal.
if strcmpi(direction,'normal')
    set(gca,'ydir','normal')
end

