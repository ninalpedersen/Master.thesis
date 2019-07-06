%Der skal lave et HSI billede til test af ICA.
%Der er tre komponenter blank, imada og newtec

%Der laves et billede blank

blank = ones(200,400);
%imwrite(blank,"blank.png")

%imada logo indlæses
circle=mean(imread("circle.png"),3);
%imada logoet gøres binært
circle(circle<200)=0;
circle(circle>=200)=1;

%newtec logo indlæses
star=mean(imread("star.png"),3);
%newtec logoet gøres binært
star(star<0.5)=0;
star(star>=0.5)=1;

x = 1:216;
circle_spektrum=sin(2*pi*7/216*x)+1;
star_spektrum=sin(2*pi*3/216*x)+1;
blank_spektrum=sin(pi/216*x)+1;

%Billeder laves nu om til HSI

temp=ones(80000,216);
temp([1:80000](circle(:)==0),:)=repmat(imada_spektrum,numel([1:80000](imada(:)==0)),1);
temp([1:80000](circle(:)==1),:)=repmat(blank_spektrum,numel([1:80000](imada(:)==1)),1);
circle=temp;

temp([1:80000](star(:)==0);:)=repmat(newtec_spektrum,numel([1:80000](newtec(:)==0)),1);
temp([1:80000](star(:)==1),:)=repmat(blank_spektrum,numel([1:80000](newtec(:)==1)),1);
star=temp;

%Der laves et tilfældigt blandingsforhold.

blandet=repmat(rand(80000,1),1,216).*star+repmat(rand(80000,1),1,216).*circle;
%blandet numeres
blandet=blandet./((blandet*ones(216,1))*ones(1,216));

%billedet af et æble indlæses
rickmorty=mean(imread("rickmorty.png"),3);
rickmorty=rickmorty/255;

blandet=blandet.*(rickmorty'*ones(1,216))/216;

%billede=GrayHyperImage(blandet);
%imshow(billede)

%save -ascii syntetisk_aeble.csv blandet


%Lad os prøve at se på data
%x=[1:216];
%plot(x,(newtec_spektrum-imada_spektrum),"b")

%vi numeret æblet væk
blandet=blandet./((blandet*ones(216,1))*ones(1,216));

billede=GrayHyperImage(blandet(:,86));
imshow(billede)



