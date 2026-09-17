function [licImage] = LIC(U, V, texture, windo, numSteps, stepSize)

[NY, NX] = size(U);
licImage = zeros(NY, NX);

% Bilinear interp
function val = interpField(F, x, y)

    ix = floor(x);
    iy = floor(y);

    fx = x - ix;
    fy = y - iy;

    ix = max(1, min(ix, NX-1));
    iy = max(1, min(NY-1, iy));

    v11 = F(iy,   ix);
    v21 = F(iy,   ix+1);
    v12 = F(iy+1, ix);
    v22 = F(iy+1, ix+1);

    val = (1-fx)*(1-fy)*v11 + ...
          fx*(1-fy)*v21 + ...
          (1-fx)*fy*v12 + ...
          fx*fy*v22;
end

% LIC LOOP
for i = 1:NX
    for j = 1:NY

        x0 = i;
        y0 = j;

        sumVal = texture(j,i);
        weightSum = windo(1);

        % Forward integration
        x = x0;
        y = y0;

        prevVec = [interpField(U,x,y), interpField(V,x,y)];

        for k = 1:numSteps

            vcurr = [interpField(U,x,y), interpField(V,x,y)];

            % normalize (important for stability)
            n = norm(vcurr);
            if n > 0
                vcurr = vcurr / n;
            end

            % sign consistency
            if dot(vcurr, prevVec) < 0
                vcurr = -vcurr;
            end

            prevVec = vcurr;

            % step
            x = x + vcurr(1)*stepSize;
            y = y + vcurr(2)*stepSize;

            %break if leave image domain
            if x < 1 || x > NX || y < 1 || y > NY
                break;
            end

            ix = floor(x);
            iy = floor(y);

            fx = x - ix;
            fy = y - iy;

            if ix >= 1 && ix < NX && iy >= 1 && iy < NY

                w = windo(k);

                value = (1-fx)*(1-fy)*texture(iy,ix) + ...
                    fx*(1-fy)*texture(iy,ix+1) + ...
                    (1-fx)*fy*texture(iy+1,ix) + ...
                    fx*fy*texture(iy+1,ix+1);

                if value >= 0.05
                    sumVal = sumVal + value*w;
                    weightSum = weightSum + w;
                end
            end
        end


        % Backward integration
        x = x0;
        y = y0;

        prevVec = [interpField(U,x,y), interpField(V,x,y)];

        for k = 1:numSteps

            vcurr = [interpField(U,x,y), interpField(V,x,y)];

            n = norm(vcurr);
            if n > 0
                vcurr = vcurr / n;
            end

            if dot(vcurr, prevVec) < 0
                vcurr = -vcurr;
            end

            prevVec = vcurr;

            x = x - vcurr(1)*stepSize;
            y = y - vcurr(2)*stepSize;

            if x < 1 || x > NX || y < 1 || y > NY
                break;
            end

            ix = floor(x);
            iy = floor(y);

            fx = x - ix;
            fy = y - iy;

            if ix >= 1 && ix < NX && iy >= 1 && iy < NY

                w = windo(k);

                value = (1-fx)*(1-fy)*texture(iy,ix) + ...
                    fx*(1-fy)*texture(iy,ix+1) + ...
                    (1-fx)*fy*texture(iy+1,ix) + ...
                    fx*fy*texture(iy+1,ix+1);

                if value >= 0.05
                    sumVal = sumVal + value*w;
                    weightSum = weightSum + w;
                end
            end
        end

        licImage(j,i) = sumVal / weightSum;

    end
end
end