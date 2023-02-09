function rv = image_coordinates(filename)
    img = imread(filename);
    % code from matlab answer forums
    % Extract the individual red, green, and blue color channels.
    redChannel = img(:, :, 1);
    greenChannel = img(:, :, 2);
    blueChannel = img(:, :, 3);
    % Find where each color plane is your exact desired color value.
    mask = redChannel == 255 & greenChannel == 0;
    [rows, columns] = find(mask);
    rv(1,:) = columns; % set up the coordinate matrix
    rv(2,:) = -rows;
%     rv = smoothdata(rv, 2, 'loess'); % smooth the sampled points. this could introduce some error
%     rv(1,:) = movmean(rv(1,:),3,2); % smooth the sampled points. this could introduce some error
%     rv(2,:) = movmean(rv(2,:),3,2);
%     rv(1,end) = max(rv(1,end-1:end)); % ensure no indent at the tip
    rv(2,:) = rv(2,:) - min(rv(2,:)); % center at (0,0)
    rv(1,:) = rv(1,:) - min(rv(1,:));
    rv = rv / max(rv(2,:)); % normalize to have height 1
end