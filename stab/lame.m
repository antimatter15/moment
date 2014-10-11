close all
clear all

filename = 'pushin.mp4';
hVideoSrc = vision.VideoFileReader(filename, 'ImageColorSpace', 'Intensity');

rgbVidplayer = vision.VideoFileReader(filename, 'ImageColorSpace', 'RGB');

% hVPlayer = vision.VideoPlayer; % Create video viewer

% poplar as a tree
% elegans as a worm
% talonted as a dinosawr

% Process all frames in the video
movMean = step(hVideoSrc);
firstrgbframe = step(rgbVidplayer);


imgB = movMean;
imgBp = imgB;
mask = ones(size(imgB,1), size(imgB,2));
correctedMean = imgBp;
ii = 2;
Hcumulative = eye(3);

numframes = 40;

hashtagB = imgB;
hashtagBp = imgBp;

frames = firstrgbframe;

kabo = mask;

while ~isDone(hVideoSrc) && ii <= numframes
    % Read in new frame
    imgA = hashtagB; % z^-1
    imgAp = hashtagBp; % z^-1
    Hcumulative = eye(3);
    
    imgB = step(hVideoSrc);
    
    rgbframe = step(rgbVidplayer);
    
    movMean = movMean + imgB;

    % Estimate transform from frame A to frame B, and fit as an s-R-t
    H = cvexEstStabilizationTform(imgA,imgB);
    HsRt = H;%cvexTformToSRT(H);
    Hcumulative = HsRt * Hcumulative;
    imgBp = imwarp(imgB,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));

    warpedrgb = imwarp(rgbframe,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));
    
    mask = mask&(imwarp(kabo,affine2d(Hcumulative),'OutputView',imref2d(size(imgB))));
    
    frames = cat(4,frames,warpedrgb);

    % Display as color composite with last corrected frame
    % step(hVPlayer, imfuse(imgAp,imgBp,'ColorChannels','red-cyan'));
    correctedMean = correctedMean + imgBp;
    
    ii = ii+1
end

little = mostrect(mask);

writerObj = VideoWriter('stabby2');
open(writerObj);

squeezy = []

little

fileID = fopen('nine.bin','w');

num=0;


for ii=1:numframes
    cropped = imresize(imcrop(frames(:,:,:,ii),little),3/18.11);
    cropped = cropped.*(cropped>0).*(cropped<1);
%     squeezy(:,:,:,ii) = cropped;
    for row=1:size(cropped,1)
        vals = uint8(255*cropped(row,:,:));
        for pixel=1:size(vals,2)
            fwrite(fileID, vals(1,pixel,:));
%             vals(1,pixel,:)
            num=num+1;
        end
    end

    whos cropped

    writeVideo(writerObj,cropped);
end

fclose(fileID);
    whos squeezy
    

    'numel squeezy'
    numel(squeezy)
    


close(writerObj)

correctedMean = correctedMean/(ii-2);
movMean = movMean/(ii-2);


'saved to file'


% Here you call the release method on the objects to close any open files
% and release memory.
release(hVideoSrc);
% release(hVPlayer);

figure; imshowpair(movMean, correctedMean, 'montage');
title(['Raw input mean', repmat(' ',[1 50]), 'Corrected sequence mean']);