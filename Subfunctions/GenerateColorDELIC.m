function colorImage = GenerateColorDELIC(OpAx, licImage, cmap)
    % Generate a color DELIC image
    colorImage = zeros(size(OpAx, 1), size(OpAx, 2), 3);  % Initialize the color image
    
    for i = 1:size(OpAx, 2)
        for j = 1:size(OpAx, 1)
            angle = abs(OpAx(j, i));  % Calculate the angle from the optical axis
            hue = angle / 180;  % Normalize the angle to [0, 1]

            % Interpolate colormap based on the normalized angle
            colorRGB = interp1(linspace(0, 1, size(cmap, 1)), cmap, hue);

            % Apply the LIC image intensity to the color channels
            colorImage(j, i, 1) = licImage(j, i) * colorRGB(1);
            colorImage(j, i, 2) = licImage(j, i) * colorRGB(2);
            colorImage(j, i, 3) = licImage(j, i) * colorRGB(3);
        end
    end
end