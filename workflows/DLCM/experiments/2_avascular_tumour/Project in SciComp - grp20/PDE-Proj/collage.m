function collImg = collage 
%#COLLIMG creates a collage of three images called 'a.png' 'b.png' and 'c.png'
%#
%# OUTPUT collImg : collage image, with individual images arranged as [a;b,c]
%#

T1 = imread('Bilder/Collage/T=1.png');
T2= imread('Bilder/Collage/T=50.png');
T3= imread('Bilder/Collage/T=100.png');
T4= imread('Bilder/Collage/T=150.png');
T5= imread('Bilder/Collage/T=201.png');

PT1 = imread('Bilder/Collage/PT=1.png');
PT2= imread('Bilder/Collage/PT=50.png');
PT3= imread('Bilder/Collage/PT=100.png');
PT4= imread('Bilder/Collage/PT=150.png');
PT5= imread('Bilder/Collage/PT=201.png');


T1 = imcrop(T1,[183 60 535 530]);
T2 = imcrop(T2,[183 60 535 530]);
T3 = imcrop(T3,[183 60 535 530]);
T4 = imcrop(T4,[183 60 535 530]);
T5 = imcrop(T5,[183 60 535 530]);

PT1 = imcrop(PT1,[140 60 530 525]);
PT2 = imcrop(PT2,[140 60 530 525]);
PT3 = imcrop(PT3,[140 60 530 525]);
PT4 = imcrop(PT4,[140 60 530 525]);
PT5 = imcrop(PT5,[140 60 530 525]);

newImageSize = [512,512]; %# or anything else that is even

%# get the new sizes - this approach requires even image size
% newSizeA = newImageSize./[2,1];
% newSizeB = newImageSize./[2,2];
% newSizeC = newImageSize./[2,2];

newSize= newImageSize./[2,2];


%# resize the images and stick together
%# place a in the top half
%# place b in the bottom left
%# place c in the bottom right 
collImg = [imresize(T1,newSize), imresize(PT1,newSize);
    imresize(T2,newSize),imresize(PT2,newSize);
    imresize(T3,newSize),imresize(PT3,newSize);
    imresize(T4,newSize),imresize(PT4,newSize);
    imresize(T5,newSize),imresize(PT5,newSize);];

% % collImg = [imresize(T1,newSize); imresize(T1,newSize);
% %     imresize(T1,newSize),imresize(T1,newSize);
% %     imresize(T1,newSize),imresize(T1,newSize);
% %     imresize(T1,newSize),imresize(T1,newSize);
% %     imresize(T1,newSize),imresize(T1,newSize);];
% 
% collImg = [imresize(T1,newSize), imresize(T1,newSize);
%     imresize(T1,newSize),imresize(T1,newSize);];

%# display the image
figure,imshow(collImg);

filename = 'Bilder/Collage/coll.png';
imwrite(collImg,filename);

