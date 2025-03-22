function [licImage] = LIC(U,V,texture,windo,numSteps,stepSize)

[NY, NX] = size(U);
licImage = zeros(NY, NX);

%Generate linear interpolant Objects
% [xGrid, yGrid] = meshgrid(linspace(1, size(scalarField,2), size(scalarField,2)), linspace(1, size(scalarField,1), size(scalarField,1)));
% interpObj_U = griddedInterpolant(yGrid, xGrid, U,'linear','none');
% interpObj_V = griddedInterpolant(yGrid, xGrid, V,'linear','none');

for i = 1:NX %Iterate through lateral pos
    for j = 1:NY %Iterate through depth pos

        % Define a boundary clipping function. Defining it within the parfor reduces overhead
        assertV = @(val, minVal, maxVal) max(min(val, maxVal), minVal);
        u = @(x, y) U(assertV(round(x), 1, NY), assertV(round(y), 1, NX));
        v = @(x, y) V(assertV(round(x), 1, NY), assertV(round(y), 1, NX));

        % assertV = @(val, minVal, maxVal) max(min(val, maxVal), minVal);
        % u = @(x, y) interpObj_U(assertV(x, 1, maxY), assertV(y, 1, maxX));
        % v = @(x, y) interpObj_V(assertV(x, 1, maxY), assertV(y, 1, maxX));

        % Initialize streamline
        x = i;
        y = j;

        % Forward integration with Gaussian window
        sum = texture(j, i) * exp(-(0^2)/(2*1^2));
        weightSum = exp(-(0^2)/(2*1^2));

        for k = 1:numSteps
            % Update position using the RK4
            % Calculate the k values
            k1U = u(y, x);
            k1V = v(y, x);

            k2U = u(y + 0.5 * k1U, x + 0.5 * k1V);
            k2V = v(y + 0.5 * k1U, x + 0.5 * k1V);

            k3U = u(y + 0.5 * k2U, x + 0.5 * k2V);
            k3V = v(y + 0.5 * k2U, x + 0.5 * k2V);

            k4U = u(y + k3U, x + k3V);
            k4V = v(y + k3U, x + k3V);

            % Apply RK4 step
            x = x + (1/6) * (k1U + 2*k2U + 2*k3U + k4U).*stepSize;
            y = y + (1/6) * (k1V + 2*k2V + 2*k3V + k4V).*stepSize;

            % Ensure streanlines are within bounds
            if x < 1 || x > NX || y < 1 || y > NY
                break;
            end

            % 2D linear interp

            ix = floor(x);
            iy = floor(y);
            fx = x - ix;
            fy = y - iy;


            if ix >= 1 && ix < NX && iy >= 1 && iy < NY
                % Calculate Gaussian weight
                weight = windo(k);

                % Accumulate texture with Gaussian weighting
                % Accumulate texture with Gaussian weighting
                value = (1 - fx) * (1 - fy) * texture(iy, ix) + ...
                    fx * (1 - fy) * texture(iy, ix + 1) + ...
                    (1 - fx) * fy * texture(iy + 1, ix) + ...
                    fx * fy * texture(iy + 1, ix + 1);

                if value < 0.05 % Helps to reduce dark regions in DELIC image 
                    continue;
                end

                sum = sum + value * weight;
                weightSum = weightSum + weight;
            end
        end

        % Backward integration with Gaussian window
        x = i;
        y = j;
        for k = 1:numSteps
            % Update position using the RK4
            k1U = u(y, x);
            k1V = v(y, x);

            k2U = u(y - 0.5 * k1U, x - 0.5 * k1V);
            k2V = v(y - 0.5 * k1U, x - 0.5 * k1V);

            k3U = u(y - 0.5 * k2U, x - 0.5 * k2V);
            k3V = v(y - 0.5 * k2U, x - 0.5 * k2V);

            k4U = u(y - k3U, x - k3V);
            k4V = v(y - k3U, x - k3V);

            % Apply RK4 step
            x = x - (1/6) * (k1U + 2*k2U + 2*k3U + k4U).*stepSize;
            y = y - (1/6) * (k1V + 2*k2V + 2*k3V + k4V).*stepSize;

            % Ensure indices are within bounds
            if x < 1 || x > NX || y < 1 || y > NY
                break;
            end

            % Bilinear interpolation
            ix = floor(x);
            iy = floor(y);
            fx = x - ix;
            fy = y - iy;
            %
            if ix >= 1 && ix < NX && iy >= 1 && iy < NY
                % Calculate Gaussian weight
                weight = windo(k);

                % Accumulate texture with Gaussian weighting
                value = (1 - fx) * (1 - fy) * texture(iy, ix) + ...
                    fx * (1 - fy) * texture(iy, ix + 1) + ...
                    (1 - fx) * fy * texture(iy + 1, ix) + ...
                    fx * fy * texture(iy + 1, ix + 1);
                % 
                if value < 0.05 %0.1
                    continue;
                end

                sum = sum + value * weight;
                weightSum = weightSum + weight;
            end
        end

        % Store the result
        licImage(j, i) = sum / weightSum;

    end
end