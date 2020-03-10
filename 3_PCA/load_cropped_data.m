% Load and crop the videos based on the filter specified
function data = load_cropped_data(vidFrames, filter, scale)
    numFrames = size(vidFrames,4);
    data = zeros(numFrames, 2);
    for j = 1:numFrames
        X = vidFrames(:,:,:,j);
        Xg = double(rgb2gray(X));
        X_cropped = Xg .* filter;
        threshold = X_cropped > scale; % light

        [Y, X] = find(threshold);
        data(j,1) = mean(X);
        data(j,2) = mean(Y);
%         subplot(1 ,2 ,1)
%         imshow(uint8((threshold * 255))); drawnow
%         subplot(1 ,2 ,2)
%         imshow(uint8(X_cropped)); drawnow
    end
end

