function [roi, rectangleHandle] = select_intensity_roi(I)
%SELECT_INTENSITY_ROI Display an intensity image and select a rectangular ROI.
%
% roi = select_intensity_roi(I)
% [roi,rectangleHandle] = select_intensity_roi(I)
%
% Draw a rectangle, adjust it, then double-click inside to accept it.
% roi = [rowFirst rowLast colFirst colLast], inclusive 1-based bounds.
% These bounds can be passed directly as processing.ROI to
% read_sparse_intensity_correlation.
%
% Select on a FULL-RESOLUTION, UNCROPPED intensity image to obtain original
% sensor coordinates. Use the reader's ComputeCorrelation=false option to
% generate this image without allocating a correlation matrix.
%
% Pixel centers inside/on the drawn rectangle are selected. The accepted
% rectangle snaps to the outer edges of those pixels and is then locked.
% Returns roi=[] if drawing is cancelled, the figure/ROI is deleted, or the
% rectangle contains no pixel centers. The accepted figure remains open.
%
% Requires Image Processing Toolbox (drawrectangle).

    validateattributes(I, {'numeric','logical'}, ...
        {'real','2d','nonempty','nonsparse'}, mfilename, 'I');
    [nRows, nCols] = size(I);
    roi = [];
    rectangleHandle = [];

    fig = figure('Name', 'Select intensity ROI', 'NumberTitle', 'off', ...
        'WindowStyle', 'normal');
    ax = axes('Parent', fig);
    imagesc(ax, I);
    axis(ax, 'image');
    xlim(ax, [0.5 nCols + 0.5]);
    ylim(ax, [0.5 nRows + 0.5]);
    colormap(ax, 'hot');
    colorbar(ax);
    xlabel(ax, 'Column');
    ylabel(ax, 'Row');
    title(ax, 'Draw a rectangle, adjust it, then double-click inside');

    try
        rectangleHandle = drawrectangle(ax, 'Color', [0 1 0], ...
            'Rotatable', false, 'DrawingArea', [0.5 0.5 nCols nRows]);
        if isempty(rectangleHandle) || ~isvalid(rectangleHandle) || ...
                isempty(rectangleHandle.Position)
            return;
        end
        wait(rectangleHandle);
    catch exception
        % Closing the window while drawing may invalidate the graphics.
        if ~isgraphics(fig) || ...
                (~isempty(rectangleHandle) && ~isvalid(rectangleHandle))
            rectangleHandle = [];
            return;
        end
        rethrow(exception);
    end
    if ~isgraphics(fig) || isempty(rectangleHandle) || ~isvalid(rectangleHandle)
        rectangleHandle = [];
        return;
    end

    position = rectangleHandle.Position; % [x y width height] in image coords
    if isempty(position) || any(position(3:4) <= 0), return; end
    colFirst = max(1, ceil(position(1)));
    colLast = min(nCols, floor(position(1) + position(3)));
    rowFirst = max(1, ceil(position(2)));
    rowLast = min(nRows, floor(position(2) + position(4)));
    if rowLast < rowFirst || colLast < colFirst
        title(ax, 'No pixel centers selected; call select_intensity_roi again');
        return;
    end
    roi = [rowFirst rowLast colFirst colLast];

    % Draw exact pixel boundaries so the confirmed selection is unambiguous.
    rectangleHandle.Position = [colFirst - 0.5, rowFirst - 0.5, ...
        colLast - colFirst + 1, rowLast - rowFirst + 1];
    rectangleHandle.InteractionsAllowed = 'none';
    title(ax, sprintf('Rows %d:%d, columns %d:%d', ...
        rowFirst, rowLast, colFirst, colLast));
end
