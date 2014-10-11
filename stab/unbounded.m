close all
clear all

filename = 'pushin.mp4';
hVideoSrc = vision.VideoFileReader(filename, 'ImageColorSpace', 'Intensity');

hVPlayer = vision.VideoPlayer; % Create video viewer

% Process all frames in the video
movMean = step(hVideoSrc);
imgB = movMean;
imgBp = imgB;
mask = ones(size(imgB,1), size(imgB,2));
correctedMean = imgBp;
ii = 2;
Hcumulative = eye(3);

numframes = 5 * 30;

frames = imgB;

kabo = mask;

while ~isDone(hVideoSrc) && ii <= numframes
    % Read in new frame
    imgA = imgB; % z^-1
    imgAp = imgBp; % z^-1
    imgB = step(hVideoSrc);
    movMean = movMean + imgB;

    % Estimate transform from frame A to frame B, and fit as an s-R-t
    H = cvexEstStabilizationTform(imgA,imgB);
    HsRt = H;%cvexTformToSRT(H);
    Hcumulative = HsRt * Hcumulative;
    imgBp = imwarp(imgB,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));

    mask = mask&(imwarp(kabo,affine2d(Hcumulative),'OutputView',imref2d(size(imgB))));
    frames = cat(4,frames,imgBp);

    % Display as color composite with last corrected frame
    % step(hVPlayer, imfuse(imgAp,imgBp,'ColorChannels','red-cyan'));
    correctedMean = correctedMean + imgBp;
    
    ii = ii+1
end

little = mostrect(mask);

writerObj = VideoWriter('stabby');
open(writerObj);
for ii=1:numframes
    writeVideo(writerObj,imcrop(frames(:,:,:,ii),little));
end
close(writerObj)

correctedMean = correctedMean/(ii-2);
movMean = movMean/(ii-2);



% Here you call the release method on the objects to close any open files
% and release memory.
release(hVideoSrc);
release(hVPlayer);

figure; imshowpair(movMean, correctedMean, 'montage');
title(['Raw input mean', repmat(' ',[1 50]), 'Corrected sequence mean']);