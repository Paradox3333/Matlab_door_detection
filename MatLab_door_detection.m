%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        Matlab Door Detection                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%%% Read in a test image with a door %%%%%

% First test image 
image = imread('./01 - R2441 - i.JPG');

% Second test image
%image = imread('./01 - R2442 - i.JPG');

[M, N, C] = size(image);
figure('Name', '01 - R2441 - i.JPG'), imshow(image,'Border','tight');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1: Image Preprocessing - Contrast Adjustment, Noise Reduction      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Conversion to Gray Value Image %%%%
imageGray = rgb2gray(image);
figure('Name','Gray image.');
imshow(imageGray,'Border','tight');

%%%% Adjust Contrast %%%%
imageContrastAdjusted = imadjust(imageGray);
figure('Name', 'Contrast-adjusted image.');
imshow(imageContrastAdjusted, 'Border', 'tight');

%%%% Noise reduction via binomial low-pass filtering %%%% 
lowPassBionomialFilter = [1,4,6,4,1]*1/16;
% 1 2D filtering operation = 2 1D filtering operations = separability.
imageBinomialFiltered = imfilter(imfilter(imageContrastAdjusted, lowPassBionomialFilter,'replicate'), lowPassBionomialFilter','replicate');
figure('Name', 'Bionomial-filtered image.');
imshow(imageBinomialFiltered, 'Border','tight');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 2: Feature Extraction - Contours                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Optimized Sobel x and y %%%%
optimizedSobel = [-3,0,3; -10,0,10; -3,0,3]/32;

gradX = imfilter(double(imageBinomialFiltered), optimizedSobel, 'replicate');
gradY = imfilter(double(imageBinomialFiltered), optimizedSobel', 'replicate');

figure('Name', 'Gradient in x-direction.');
imagesc(gradX);
colormap hsv;
axis off;
colorbar;

figure('Name', 'Gradient in y-direction.');
imagesc(gradY);
colormap hsv;
axis off;
colorbar;

gradAbsX = sqrt(double(gradX.*gradX));
gradAbsY = sqrt(double(gradY.*gradY));

figure('Name', 'Gradient magnitude in x-direction, in log scale.');
imagesc(log(double(gradAbsX+1)));
colormap parula;
colorbar;
axis off;

figure('Name', 'Gradient magnitude in y-direction, in log scale.');
imagesc(log(double(gradAbsY+1)));
colormap parula;
colorbar;
axis off;

%%%% Edge Strength %%%%

% Negative x-Gradients (bright to dark)  %
positionMapNegGradX = (gradX < 0);
negGradXdupe = gradX .* positionMapNegGradX;
negGradXdupe(~positionMapNegGradX) = nan;

% Define edge strength threshold %
negGradXmean = mean(negGradXdupe,"all",'omitnan');

% Segment strong negative x-Gradients %
thresholdNegGradX = (gradX < 2.5*negGradXmean);

figure('Name', 'Strong light-dark transitions (negative gradient) in x-direction, as binary image.');
imagesc(1-thresholdNegGradX);
colormap gray;
colorbar;
axis off;

% Noise reduction %
noiseReducedNegGradX = medfilt2(thresholdNegGradX);
figure('Name', 'Noise-reduced strong light-dark transitions (negative gradient) in x-direction.');
imagesc(1-noiseReducedNegGradX);
colormap gray;
colorbar;
axis off;

% Dilate for 1 pixel %
structuralElement1 = strel('square',4);
structuralElement2 = strel('square',4);

morphedNegGradX = imclose(imdilate(noiseReducedNegGradX, structuralElement1), structuralElement2);
figure('Name', 'Continuous transitions + small size transitions removed, for negative gradient in x-direction.');
imagesc(1-morphedNegGradX);
colormap gray;
colorbar;
axis off;

% Positive x-Gradients (dark to bright) %
positionMapPosGradX = (gradX > 0);
posGradXdupe = gradX .* positionMapPosGradX;
posGradXdupe(~positionMapPosGradX) = nan;
posGradXmean = mean(posGradXdupe,"all",'omitnan');

% Segment strong positive x-Gradients %
thresholdPosGradX = (gradX > 2.5*posGradXmean);

figure('Name', 'Strong dark-light transitions (positive gradient) in x-direction, as binary image.');
imagesc(1-thresholdPosGradX);
colormap gray;
colorbar;
axis off;

% Noise reduction %
noiseReducedPosGradX = medfilt2(thresholdPosGradX);
figure('Name', 'Noise-reduced strong dark-light transitions (positive gradient) in x-direction.');
imagesc(1-noiseReducedPosGradX);
colormap gray;
colorbar;
axis off;

% Dilate for 1 pixel %
morphedPosGradX = imclose(imdilate(noiseReducedPosGradX, structuralElement1), structuralElement2);
figure('Name', 'Continuous transitions + small size transitions removed, for positive gradient in x-direction');
imagesc(1-morphedPosGradX);
colormap gray;
colorbar;
axis off;

% Negative y-Gradients (bright to dark)  %
positionMapNegGradY = (gradY < 0);
negGradYdupe = gradY .* positionMapNegGradY;
negGradYdupe(~positionMapNegGradY) = nan;

% Define y-edge strength threshold %
negGradYmean = mean(negGradYdupe,"all",'omitnan');

% Segment strong negative y-Gradients %
thresholdNegGradY = (gradY < 3*negGradYmean);

figure('Name', 'Strong light-dark transitions (negative gradient) in y-direction, as binary image.');
imagesc(1-thresholdNegGradY);
colormap gray;
colorbar;
axis off;

% Noise reduction %
noiseReducedNegGradY = medfilt2(thresholdNegGradY);
figure('Name', 'Noise-reduced strong light-dark transitions (negative gradient) in y-direction.');
imagesc(1-noiseReducedNegGradY);
colormap gray;
colorbar;
axis off;

% Dilate for 1 pixel %
morphedNegGradY = imclose(imdilate(noiseReducedNegGradY, structuralElement1), structuralElement2);
figure('Name', 'Continuous transitions + small size transitions removed, for negative gradient in y-direction.');
imagesc(1-morphedNegGradY);
colormap gray;
colorbar;
axis off;

% Positive y-Gradients (dark to bright) %
positionMapPosGradY = (gradY > 0);
posGradYdupe = gradY .* positionMapPosGradY;
posGradYdupe(~positionMapPosGradY) = nan;
posGradYmean = mean(posGradYdupe,"all",'omitnan');

% Segment strong positive y-Gradients %
thresholdPosGradY = (gradY > 2.6*posGradYmean);

figure('Name', 'Strong dark-light transitions (positive gradient) in y-direction, as binary image.');
imagesc(1-thresholdPosGradY);
colormap gray;
colorbar;
axis off;

% Noise reduction %
noiseReducedPosGradY = medfilt2(thresholdPosGradY);
figure('Name', 'Noise-reduced strong dark-light transitions (positive gradient) in y-direction.');
imagesc(1-noiseReducedPosGradY);
colormap gray;
colorbar;
axis off;

% Dilate for 1 pixel %
structuralElement3 = strel('square',1);
structuralElement4 = strel('square',1);

% Dilate for 1 pixel %
morphedPosGradY = imclose(imdilate(noiseReducedPosGradY, structuralElement3), structuralElement4);
figure('Name', 'Continuous transitions + small size transitions removed, for positive gradient in y-direction.');
imagesc(1-morphedPosGradY);
colormap gray;
colorbar;
axis off;

%load('Solutions_Task_2.mat')

%%%% Visualize Results %%%%
% figure(12),imagesc(img_sobel_x); colorbar; colormap hsv; axis off;
% figure(13),imagesc(log(1+img_edge_strength_x)); colorbar; axis off;
% figure(14),imagesc(img_neg_edge_x); colorbar; colormap gray; axis off;
% figure(15),imagesc(img_pos_edge_x); colorbar; colormap gray; axis off;
% figure(16),imagesc(img_neg_edge_y); colorbar; colormap gray; axis off;
% figure(17),imagesc(img_pos_edge_y); colorbar; colormap gray; axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 3: Segmentation of Door Gap Contour                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Vertical Edge candidates for %
% negPosGradXcombined = (img_pos_edge_x & img_neg_edge_x);
negPosGradXcombined = (morphedPosGradX & morphedNegGradX);
negPosGradXcombined = medfilt2(negPosGradXcombined);

% Design a symmetric filter %
negPosGradXcombinedFiltered = imfilter(negPosGradXcombined, fspecial('average'));

% Combine knowledge on gray value range and bright-to-dark gradients 
negPosGradXcombinedFiltered = medfilt2(negPosGradXcombinedFiltered);

% Fatten lines to improve detection 
xFinal = imclose(imdilate(negPosGradXcombinedFiltered, structuralElement1), structuralElement2);

figure('Name', 'Negative and positive gradients-X combined and filtered.');
imagesc(1 - xFinal);
colormap gray;
colorbar;
axis off;

% Extract Horizontal Edge candidates 
negPosGradYcombined = (morphedPosGradY | morphedNegGradY);
negPosGradYcombined = medfilt2(negPosGradYcombined);

% Combine knowledge on absolute gray value and bright-to-dark-gradient
negPosGradYcombinedFiltered = imfilter(negPosGradYcombined, fspecial('average')');
negPosGradYcombinedFiltered = medfilt2(negPosGradYcombinedFiltered);

% Fatten lines to improve detection 
yFinal = imclose(imdilate(negPosGradYcombinedFiltered, structuralElement3), structuralElement4);

figure('Name', 'Negative and positive gradients-Y combined and filtered.');
imagesc(1 - yFinal);
colormap gray;
colorbar;
axis off;

% load('Solutions_Task_3.mat')

%%%% Visualize Results %%%%
% figure(8),imagesc(max(max(img_sym_x))-abs(img_sym_x));colorbar; colormap gray;axis off;
% figure(9),imagesc(1-candidates_x);colorbar; colormap gray;axis off;
% figure(10),imagesc(1-candidates_y);colorbar; colormap gray;axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 4: Measuring Lines of the Door Gap                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Hough space for vertical lines
[H1, T1, R1] = hough(xFinal);

% hPeaksX3lines = houghpeaks(H1, 3); Matlab is bugged! We can't update the value (2 or 3) as it forever gives the output of the first run of the function. 

hPeaksX = houghpeaks(H1, 2);

figure('Name', "Hough space for vertical lines. ");
imagesc(H1, 'XData', T1, 'YData', R1);
colormap gray;
hold on;
plot(T1(hPeaksX(:,2)), R1(hPeaksX(:,1)), 's', 'color', 'white');

%%%% Find vertical Hough Lines %%%%
verticalHoughLines = houghlines(xFinal, T1, R1, hPeaksX, 'FillGap', 1000, 'MinLength', 5);
figure('Name', 'Vertical Hough lines.');
imshow(imageGray);
hold on;

for k = 1:length(verticalHoughLines)
    xy = [verticalHoughLines(k).point1; verticalHoughLines(k).point2]; 
    plot(xy(:,1), xy(:,2), 'LineWidth', 1, 'Color','green');
end

% Transform two best candidates to Hessian Normal form
for k = 1:2
    hessianLinesX(:,k) = [cos(pi/180 * verticalHoughLines(k).theta); sin(pi/180 * verticalHoughLines(k).theta); -verticalHoughLines(k).rho];
end

% Compute Hough space for horizontal lines
[H2, T2, R2] = hough(yFinal);
hPeaksY = houghpeaks(H2, 2);

figure('Name', "Hough space for horizontal lines. ");
imagesc(H2, 'XData', T2, 'YData', R2);
colormap gray;
hold on;
plot(T1(hPeaksY(:,2)), R1(hPeaksY(:,1)), 's', 'color', 'white');

%%%% Find horizontal Hough Lines %%%%
horizontalHoughLines = houghlines(yFinal, T2, R2, hPeaksY, 'FillGap', 1000, 'MinLength', 5);
figure('Name', 'Horizontal Hough lines.');
imshow(imageGray);
hold on;

for k = 1:length(horizontalHoughLines)
    xy = [horizontalHoughLines(k).point1; horizontalHoughLines(k).point2]; 
    plot(xy(:,1), xy(:,2), 'LineWidth', 1, 'Color','green');
end

% Transform two best candidates to Hessian Normal form
for k = 1:2
    hessianLinesY(:,k) = [cos(pi/180 * horizontalHoughLines(k).theta); sin(pi/180 * horizontalHoughLines(k).theta); -horizontalHoughLines(k).rho];
end

% Plot lines
figure('Name', 'Vertical and horizontal Hough lines of interest.'), imshow(image,'Border','tight'); hold on;

for k=1:2
    plot([round(-(hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k)); round(-(M*hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k));], [1,M], 'LineWidth',1,'Color','green');
end

% plot a line through 2 points, through the function 
% x*cos(theta) + y*sin(theta) = rho.

for k=1:2
    plot([1,N], [round(-(hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k)); round(-(N*hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k));], 'LineWidth',1,'Color','green');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 5: Measure and Classify Corner Points                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find left and right door gap lines %
figure('Name', 'Distinction between left and right Hough lines.'), imshow(image,'Border','tight'); hold on;

% Right
for k=1
    plot([round(-(hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k)); round(-(M*hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k));], [1,M], 'LineWidth',1,'Color','green');
end

% Left
for k=2
    plot([round(-(hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k)); round(-(M*hessianLinesX(2,k)+hessianLinesX(3,k))/hessianLinesX(1,k));], [1,M], 'LineWidth',1,'Color','red');
end

hold off;

% Find up and down door gap lines %
figure('Name', 'Distinction between top and bottom Hough lines.'), imshow(image,'Border','tight'); hold on;

% Bottom
for k=1
    plot([1,N], [round(-(hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k)); round(-(N*hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k));], 'LineWidth',1,'Color','green');
end

% Top
for k=2
    plot([1,N], [round(-(hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k)); round(-(N*hessianLinesY(1,k)+hessianLinesY(3,k))/hessianLinesY(2,k));], 'LineWidth',1,'Color','red');
end

hold off;

% Find corner points %
for k = 1:M
    leftLinePoints(1,k) = round(-(k*hessianLinesX(2,2)+hessianLinesX(3,2))/hessianLinesX(1,2));
    leftLinePoints(2,k) = k;

    rightLinePoints(1,k) = round(-(k*hessianLinesX(2,1)+hessianLinesX(3,1))/hessianLinesX(1,1));
    rightLinePoints(2,k) = k;
end

for k = 1:N
    bottomLinePoints(2,k) = round(-(k*hessianLinesY(1,1)+hessianLinesY(3,1))/hessianLinesY(2,1));
    bottomLinePoints(1,k) = k;

    topLinePoints(2,k) = round(-(k*hessianLinesY(1,2)+hessianLinesY(3,2))/hessianLinesY(2,2));
    topLinePoints(1,k) = k;
end

figure('Name', 'Line intersection points.'), imshow(image,'Border','tight'); hold on;
plot(bottomLinePoints(1,:), bottomLinePoints(2,:), 'LineWidth',1,'Color','green');
plot(topLinePoints(1,:), topLinePoints(2,:), 'LineWidth',1,'Color','green');
plot(rightLinePoints(1,:), rightLinePoints(2,:), 'LineWidth',1,'Color','green');
plot(leftLinePoints(1,:), leftLinePoints(2,:), 'LineWidth',1,'Color','green');

topLeftLinesIntersect = intersect(leftLinePoints', topLinePoints', "rows");
plot(topLeftLinesIntersect(1,1), topLeftLinesIntersect(1,2),'go', 'MarkerFaceColor','g');

topRightLinesIntersect = intersect(rightLinePoints', topLinePoints', "rows");
plot(topRightLinesIntersect(1,1), topRightLinesIntersect(1,2),'bo', 'MarkerFaceColor','b');

bottomRightLinesIntersect = intersect(rightLinePoints', bottomLinePoints', "rows");
plot(bottomRightLinesIntersect(1,1), bottomRightLinesIntersect(1,2),'ko', 'MarkerFaceColor','k');

bottomLeftLinesIntersect = intersect(bottomLinePoints', leftLinePoints', "rows");
plot(bottomLeftLinesIntersect(1,1), bottomLeftLinesIntersect(1,2),'ro', 'MarkerFaceColor','r');