
%function to find the background of an image, I, using a 4th order
%polynomial fit to the background pixels
%input: I, the grayscale image to find the background of
%output: B, the background of I
%method: find 'objects' in I, mask them from the image, paint the remaining
%area using the inpaint_nans function

function [BWfinal,SS] = imagebackground_poly4(I)

% find high frequency features in image using sobel filter
[junk threshold] = edge(I, 'sobel');
fudgeFactor = 1; %was 0.4 for RBCs
BWs = edge(I,'sobel', threshold * fudgeFactor);

% dilate the image mask, values defined works for 120X phase images acquired in our lab
se90 = strel('line', 6, 90);
se0 = strel('line', 6, 0);
BWsdil = imdilate(BWs, [se90 se0]);

% fill gaps
BWdfill = imfill(BWsdil, 'holes');

% smooth image
seD = strel('diamond',5);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000); 

% dilate mask further
se90 = strel('line', 30, 90);
se0 = strel('line', 30, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);

IList = I(~BWfinal);

% initiate the grid in image background region
sz = size(I);
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
XList = X(~BWfinal);
YList = Y(~BWfinal);

% fit the grid surface using polyfiting
if sum(~isnan(BWfinal)) ~=0
    CFit = polyfitn([XList,YList], IList, 4);
    B =(reshape(polyvaln(CFit, [X(:), Y(:)]), sz(1), sz(2)));
else
    B = I;
end
SS=B-I;
