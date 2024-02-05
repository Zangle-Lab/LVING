%% LIVING supplementary figures
clear all; clc;

%% figure S2: Verification of QPI system baseline performance using bead data

load('T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Mar2022\WS_S2_rev1.mat');

%% figure S2a : QPI image of polystyrene beads

load('M:\Data\Soorya\stage_test\Test_accuracy_precision_melted_beads\Trigger_stageloop\Image_620.mat');
[M,Im] = imagebackground_poly4_bead(Phase);
figure(1); imagesc(Im(1:170,246:415)); axis off image;
barsize=10; pxlsize=0.76/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);

%% figure S2b : actual optical volume of beads vs QPI computed OV

figure(2);
plot(OV_calc,'*b') 
hold on;
plot(OV_ori,'k');

%% figure S2c : OV accuracy box plot

figure(3); boxplot(accuracy);

%% figure S2d : OV precision box plot

figure(4); boxplot(precision);

%% figure S3a: Velocity map of fixed RPE cell moving on stage computed by QPV

load('T:\Data\Soorya\FixedCell120XImaging_April2021\Trial1\Cell_12_pos_121.mat');
[B,Abkg1]=imagebackground_poly4_double(Phase);
load('T:\Data\Soorya\FixedCell120XImaging_April2021\Trial1\Cell_12_pos_131.mat');
[B,Abkg2]=imagebackground_poly4_double(Phase);
gs=15; xcg=4; sz=size(Abkg1);
XCs=SSD_corr_rev4(Abkg2,Abkg1,gs); 
XCsp=zeros((2*gs)+1,(2*gs)+1,sz(1),sz(2),'single');
for uu=1:sz(1)-xcg+1    %spatial averaging of SSD
    for vv=1:sz(2)-xcg+1
        XCsp(:,:,uu,vv)=mean(mean(XCs(:,:,uu:uu+xcg-1,vv:vv+xcg-1),3),4);
    end
end
clear XCs;
dX1=zeros(sz(1)+1);
dY1=dX1;
for ii =1:sz(1)-3
   for jj =1:sz(2)-3
      [xm,ym]=findvalley_v3_rev4(XCsp(:,:,ii,jj));
      dX1(ii+2,jj+2) = xm;
      dY1(ii+2,jj+2) = ym;
   end
end
clear XCsp;
medfiltsz=16; kk=1; tcg=1;
for ee=medfiltsz+1:sz(1)-medfiltsz+1
    for ff=medfiltsz+1:sz(2)-medfiltsz+1
        Ux=dX1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
        ux=median(Ux(:),'omitnan');
        Uy=dY1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
        uy=median(Uy(:),'omitnan'); 
        if abs(dX1(ee,ff))>(1.5*ux) || abs(dX1(ee,ff))<(0.7*ux)
            dX(ee,ff,kk-tcg+1)=ux;
        else
            dX(ee,ff,kk-tcg+1)=dX1(ee,ff);
        end
        if abs(dY1(ee,ff))>(1.5*uy) || abs(dY1(ee,ff))<(0.7*uy)
            dY(ee,ff,kk-tcg+1)=uy;
        else
            dY(ee,ff,kk-tcg+1)=dY1(ee,ff);
        end
    end
end
figure(3);
dX2=imresize(dX(50:300,120:370),[250/8,250/8]); 
dY2=imresize(dY(50:300,120:370),[250/8,250/8]);
[Vxo,Vyo]=meshgrid(-136:8:114,-206:8:44); %(-61:8:169,-89:8:165)
x=[-138 116];
y=[-208 46];
[M,~]=imagebackground_poly4(-Abkg1);
MB=imresize(M(50:300,120:370),0.125);
DX=MB.*dX2;
DY=MB.*dY2;
DX(31,2)=0.238*8*5;
imagesc(x,y,Abkg1(50:300,120:370)*(0.623/0.18)); colormap(temp); %(167:421,195:425)
axis off image; set(gcf,'color','w'); caxis([-0.1 1.3]);
hold on;
quiver(Vxo,Vyo,DX,DY,3,'k');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=80; ybase=30; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[90 180],'fontsize',15,'color','w','string',textdown)
textup=sprintf('5 µm');
text('units','pixels','position',[200 90],'fontsize',15,'color','k','string',textup)

%% figure S3 b : RPE fixed cell LVING growth map overlayed on QPI image

load('T:\Data\Soorya\FixedCell120XImaging_April2021\MassGenResults_rev136\Results_2hr\WS2_cell12.mat')
SGf=zeros(512,512);
SGfin = imfilter((GC), fspecial('gaussian', [50 50], 1));
SGf(5:508,5:508)=SGfin;
SM = Abkg1>0.001;
SGf = SGf.*SM;
SGf(SGf==0) = NaN;
figure(4);
imoverlay(Abkg1(50:300,120:370),SGf(50:300,120:370).*60,[-0.1, 0.1],[],parula, 0.2, gca); colormap(map);
colorbar; axis image; set(gcf,'color','w');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)

%% figure S3c: RPE fixed cell LVING vs QPI specific growth rate plot

load('T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Mar2022\WS_S3_rev1.mat')
sgrWhole(isnan(sgrWhole))=0;
h=bar([mean(sgrWhole(sgrWhole~=0)),mean(sgrFx(sgrFx~=0))]);
set(gca, 'XTickLabel', {'QPI','LIVING'});
hold on;
errorbar([1,2],[mean(sgrWhole(sgrWhole~=0)),mean(sgrFx(sgrFx~=0))],...
    [std(sgrWhole(sgrWhole~=0))/sqrt(nnz(sgrWhole)),std(sgrFx(sgrFx~=0))/sqrt(nnz(sgrFx))],...
    [std(sgrWhole(sgrWhole~=0))/sqrt(nnz(sgrWhole)),std(sgrFx(sgrFx~=0))/sqrt(nnz(sgrFx))],...
    '.','Color','b');
hold on;
ydata=[nonzeros(sgrWhole),nonzeros(mean(sgrFx,2))]; [r, c] = size(ydata); xdata = repmat(1:c, r, 1);
scatter(xdata(:), ydata(:), 'k.', 'jitter','on', 'jitterAmount', 0.05);
set(gcf,'Color','w'); ylim([-0.1 0.1]);
ylabel('Specific growth rate (/hr)');

%% figure S4: LVING mass conservation validation in live RPE cell, 
% comparison of specific growth rate between QPI and LVING

for mm=1:30
    fdir='S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\MassGenResults_rev136\Results_2hr\';
    fnameM=sprintf('WS1_cell%d.mat',mm);
    load([fdir fnameM]);
    for uu=1:120
        Bmask=Abkg_stored2(:,:,uu)>0.001;
        mass(uu)=sum(sum(Abkg_stored2(:,:,uu).*Bmask));
    end
    fitData=polyfit(Time,mass',1);
    sgrQPI(mm)=(fitData(1)/mean(mass))*60;
    [BWfinal,SS] = imagebackground_poly4(-Abkg_stored2(:,:,1));
    sgrFx(mm)=(sum(sum(GC.*BWfinal(5:508,5:508)))/mean(mass))*60;
end

%% figure S4a : Plot of dry mass vs time for RPE cell between QPI & LVING

plot(MassAm); hold on; plot(mass);
ylim([250 400]);

%% figure S4b : Comparison plot of specific growth rate between QPI & LVING

h=bar([mean(sgrQPI(sgrQPI~=0)),mean(sgrFx(sgrFx~=0))]);
set(gca, 'XTickLabel', {'QPI','LVING'});
hold on;
errorbar([1,2],[mean(sgrQPI(sgrQPI~=0)),mean(sgrFx(sgrFx~=0))],...
    [std(sgrQPI(sgrQPI~=0))/sqrt(nnz(sgrQPI)),std(sgrFx(sgrFx~=0))/sqrt(nnz(sgrFx))],...
    [std(sgrQPI(sgrQPI~=0))/sqrt(nnz(sgrQPI)),std(sgrFx(sgrFx~=0))/sqrt(nnz(sgrFx))],...
    '.','Color','b');
hold on;
ydata=[nonzeros(sgrQPI),nonzeros(sgrFx)]; [r, c] = size(ydata); xdata = repmat(1:c, r, 1);
scatter(xdata(:), ydata(:), 'r.', 'jitter','on', 'jitterAmount', 0.05);
set(gcf,'Color','w'); ylim([-0.1 0.1]);
ylabel('Specific growth rate (/hr)');


%% figure S5 : 

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
  
load('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\MassGenTrackingResults_rev136_Trial1\cell7\MassGenResults_rev13_6_trial2\WS10_cell7.mat')
SGf=zeros(512,512); SGf = imfilter((GC.*60), fspecial('gaussian', [50 50], 1));
SM = Abkg_stored2(5:508,5:508,1)>-0.01; SGf = SGf.*SM; SGf(SGf==0) = NaN;
figure(1);
imoverlay(Abkg_stored2(5:508,5:508,1),SGf,[-0.03, 0.03],[-0.005 0.035],parula, 0.2, gca); 
colormap(map); ylim([70 370]); xlim([128 428]);
GFP=imread('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\Trial1\GFP120X_7_frame_570.tif');
RFP=imread('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\Trial1\RFP120X_7_frame_570.tif');
minr=100; ming=200; maxr=500; maxg=1000;
RFPmask1=RFP>minr; RFPmask2=RFP<maxr;
RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr);
GFPmask1=GFP>ming;
GFPmask2=GFP<maxg;
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
GFP2=(GFP1-ming)/(maxg-ming);
DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
figure(2)
imshow(DImage3); hold on; 
barsize=10; pxlsize=0.23/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
ylim([70 370]); xlim([128 428]);

%% figure S6 a to c: cell cytoplasm and nucleus segmentation validation

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;

load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\Trial1\QPM120X_3_frame_3.mat')
[junk threshold] = edge(Phase, 'sobel');
fudgeFactor = 1; %was 0.4 for RBCs
BWs = edge(Phase,'sobel', threshold * fudgeFactor);
se90 = strel('line', 6, 90);
se0 = strel('line', 6, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
seD = strel('diamond',5);
%BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);
[BnC,Lnlls] = bwboundaries(BWfinal,'noholes');

GFP0=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\Trial1\mAG120X_3_frame_3.tif');
RFP0=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\Trial1\mKO2120X_3_frame_3.tif');
GFP=GFP0(:,1:671); RFP=RFP0(:,1:671); 
ming=250; maxg=1000;
image=zeros(512,671,3);
minr=single(min(RFP(:)))+350;
maxr=single(max(RFP(:)))+50; 

RFPmask1=RFP>minr;
RFPmask2=RFP<maxr;
RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
image(:,:,1)=RFP2;
Image=rgb2gray(image);
%     Image=RFP2;
BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
sz=size(BWn);
c = kmeans(BWn(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
jk=mode(CC(:));
BWsn=CC~=jk;
%     BWsn = imbinarize(BWn,0.15);
seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
BWfinaln = imerode(BWsn,seD);
BWfinaln = bwareaopen(BWfinaln, 200);
Ln1 = imclearborder(BWfinaln, 4);
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ln1), kernel, 'same');
binaryImage = blurryImage > 0.5; 
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);
% binaryImage = imerode(binaryImage, [se90 se0]);
% Lnls=imresize(binaryImage,[356,458]);
Lnls=imresize(BWfinaln,[356,458]);
BWint=zeros(512,512); BWint(131:486,55:512)=Lnls;
[Bn,Lnlls] = bwboundaries(BWint,'noholes');
[~,Abkg]=imagebackground_poly4(Phase);
figure(1); imagesc(Abkg); colormap(temp); axis off image;
hold on
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
for ee=1:length(BnC)
    boundary = BnC{ee};
    plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
xlim([150 430]); ylim([150 400]);

GFP2=imresize(GFP,[356,458]);
GFP3=zeros(512,512); GFP3(131:486,55:512)=GFP2;
mycolors1=zeros(2001,3);
mycolors1(1:2001,2)=0:0.0005:1;
figure(2)
imagesc(GFP3); caxis([ming maxg]); colormap(mycolors1); hold on; 
hold on
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
for ee=1:length(BnC)
    boundary = BnC{ee};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
axis off; xlim([150 430]); ylim([150 400]);

RFP21=imresize(RFP2,[356,458]);
RFP3=zeros(512,512); RFP3(131:486,55:512)=RFP21;
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
figure(3)
imagesc(RFP3); colormap(mycolors2); hold on; 
hold on
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
for ee=1:length(BnC)
    boundary = BnC{ee};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
axis off; xlim([150 430]); ylim([150 400]);


% figure S6 d: QPV nucleus tracking validation
% load('K:\Data\Soorya\RPEFUCCIImaging_2020\Pn7_28Feb2020\MassGenResults_rev137\MassGenResults_rev137_30min\WS1_cell1.mat');
load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\MassGenResults_rev137\MassGenResults_rev137_30min\WS1_cell16.mat')
imagesc(Abkg_stored2(:,:,1)/(0.21*0.21)); axis off; colormap(temp); caxis([-0.1 1]);
GFP0=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\Trial1\mAG120X_16_frame_16.tif');
RFP0=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\Trial1\mKO2120X_16_frame_16.tif');
GFP0=GFP0(:,1:671); RFP0=RFP0(:,1:671); 
GFP30=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\Trial1\mAG120X_16_frame_46.tif');
RFP30=imread('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\Trial1\mKO2120X_16_frame_46.tif');
GFP30=GFP30(:,1:671); RFP30=RFP30(:,1:671); 

% label at time 0 minutes as starting shape of the nucleus
    % GFP=GFP0(:,1:671); RFP=RFP0(:,1:671);
    image=zeros(512,671,3);
    % mycolors1=zeros(2001,3);
    % mycolors1(1:2001,2)=0:0.0005:1;
    % %     maxg=single(min(GFP(:))+1000);
    % ming=single(min(GFP(:))+250);
    % maxg=single(max(GFP(:)))+30; 
    % %     maxg=2000; 
    minr=single(min(RFP0(:)))+150;
    maxr=single(max(RFP0(:)))+50; 
    % %     maxr=2000; minr=single(min(RFP(:))+80);
    % %     maxr=1000; %minr=100;
    % GFPmask1=GFP>ming;
    % GFPmask2=GFP<maxg;
    % GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
    % GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
    % image(:,:,2)=GFP2;

RFPmask1=RFP0>minr;
RFPmask2=RFP0<maxr;
RFP1=(single(RFP0).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,1)=RFP2;
Image=rgb2gray(image);
%     Image=RFP2;
BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
sz=size(BWn);
c = kmeans(BWn(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
jk=mode(CC(:));
BWsn=CC~=jk;
%     BWsn = imbinarize(BWn,0.15);
seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
BWfinaln = imerode(BWsn,seD);
BWfinaln = bwareaopen(BWfinaln, 200);
Ln1 = imclearborder(BWfinaln, 4);
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ln1), kernel, 'same');
binaryImage = blurryImage > 0.5; 
se90 = strel('line', 12, 90);
se0 = strel('line', 12, 0);
%     binaryImage = imdilate(binaryImage, [se90 se0]);
%     Lnls=imresize(binaryImage,[356,458]);
Lnls=imresize(BWfinaln,[356,458]);
BWint=zeros(512,512); BWint(131:486,55:512)=Lnls;
[Bn,Lnlls] = bwboundaries(BWint,'noholes');
figure(1); imagesc(Abkg_stored2(:,:,1)); colormap(temp); axis off image;
hold on
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end

% lagrangian track points
zln=30; tt=258; ww=271;
hold on   
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=258; ww=290;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=298; ww=266;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=275; ww=261;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=278; ww=302;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)

tt=313; ww=284;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=240; ww=288;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=251; ww=284;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=259; ww=283;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=268; ww=284;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=277; ww=288;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=282; ww=294;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=278; ww=310;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=265; ww=320;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=255; ww=326;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=244; ww=330;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=273; ww=315;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)
tt=233; ww=331;
xos(1,:)=reshape(XS(tt,ww,1:zln),1,zln); yos(1,:)=reshape(YS(tt,ww,1:zln),1,zln);
plot(xos,yos,'b', 'MarkerSize', 5, 'LineWidth', 2)

% label at time 30 minutes as final shape of the nucleus
    % GFP=GFP30(:,1:671); RFP=RFP30(:,1:671);
    image=zeros(512,671,3);
    % mycolors1=zeros(2001,3);
    % mycolors1(1:2001,2)=0:0.0005:1;
    % %     maxg=single(min(GFP(:))+1000);
    % ming=single(min(GFP(:))+250);
    % maxg=single(max(GFP(:)))+30; 
    % %     maxg=2000; 
    minr=single(min(RFP30(:)))+150;
    maxr=single(max(RFP30(:)))+50; 
    % %     maxr=2000; minr=single(min(RFP(:))+80);
    % %     maxr=1000; %minr=100;
    % GFPmask1=GFP>ming;
    % GFPmask2=GFP<maxg;
    % GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
    % GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
    % image(:,:,2)=GFP2;

RFPmask1=RFP30>minr;
RFPmask2=RFP30<maxr;
RFP1=(single(RFP30).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,1)=RFP2;
Image=rgb2gray(image);
%     Image=RFP2;
BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
sz=size(BWn);
c = kmeans(BWn(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
jk=mode(CC(:));
BWsn=CC~=jk;
%     BWsn = imbinarize(BWn,0.15);
seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
BWfinaln = imerode(BWsn,seD);
BWfinaln = bwareaopen(BWfinaln, 200);
Ln1 = imclearborder(BWfinaln, 4);
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ln1), kernel, 'same');
binaryImage = blurryImage > 0.5; 
se90 = strel('line', 12, 90);
se0 = strel('line', 12, 0);
%     binaryImage = imdilate(binaryImage, [se90 se0]);
%     Lnls=imresize(binaryImage,[356,458]);
Lnls=imresize(BWfinaln,[356,458]);
BWfnl=zeros(512,512); BWfnl(131:486,55:512)=Lnls;
[Bnf,Lnllf] = bwboundaries(BWfnl,'noholes');
% figure(2); imagesc(Abkg_stored2(:,:,30)); colormap(temp); axis off image;
hold on
for ee=1:length(Bnf)
    boundary = Bnf{ee};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end

%% figure S7 : 
% S7.a: MCF7 cell cycle phase growth rate

load('F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\MassGenResults_rev136\SegmentedWSs\NCWsgr_rev1.mat')

sgrN0_all = sgrN0; sgrC0_all=sgrC0; sgrW0_all=sgrW0;
sgrN0_all(isinf(abs(sgrN0_all))|isnan(sgrN0_all))=0; sgrC0_all(isinf(abs(sgrC0_all))|isnan(sgrC0_all))=0; sgrW0_all(isinf(abs(sgrW0_all))|isnan(sgrW0_all))=0;

G1WS=(CellCycle==1); G1SWS=(CellCycle==2); SG2WS=(CellCycle==3); 
SRNcG1_all=G1WS.*sgrN0_all; SRNcG1S_all=G1SWS.*sgrN0_all; SRNcSG2_all=SG2WS.*sgrN0_all; 
SRCyG1_all=G1WS.*sgrC0_all; SRCyG1S_all=G1SWS.*sgrC0_all; SRCySG2_all=SG2WS.*sgrC0_all; 
SRWG1_all=G1WS.*sgrW0_all; SRWG1S_all=G1SWS.*sgrW0_all; SRWSG2_all=SG2WS.*sgrW0_all; 

for cc=1:length(CellCycle)
    SRNcG1(cc,1)=mean(nonzeros(SRNcG1_all(cc,:))); 
    SRNcG1S(cc,1)=mean(nonzeros(SRNcG1S_all(cc,:)));
    SRNcS(cc,1)=mean(nonzeros(SRNcSG2_all(cc,:)));  

    SRCyG1(cc,1)=mean(nonzeros(SRCyG1_all(cc,:))); 
    SRCyG1S(cc,1)=mean(nonzeros(SRCyG1S_all(cc,:)));
    SRCyS(cc,1)=mean(nonzeros(SRCySG2_all(cc,:))); 
    
    SRWG1(cc,1)=mean(nonzeros(SRWG1_all(cc,:)));
    SRWG1S(cc,1)=mean(nonzeros(SRWG1S_all(cc,:)));
    SRWS(cc,1)=mean(nonzeros(SRWSG2_all(cc,:))); 
   
end
SRNcG1(isnan(SRNcG1))=0; SRNcG1S(isnan(SRNcG1S))=0; SRNcS(isnan(SRNcS))=0; 
SRCyG1(isnan(SRCyG1))=0; SRCyG1S(isnan(SRCyG1S))=0; SRCyS(isnan(SRCyS))=0; 
SRWG1(isnan(SRWG1))=0; SRWG1S(isnan(SRWG1S))=0; SRWS(isnan(SRWS))=0; 

NdifG1mn=median(nonzeros(SRNcG1)); NdifG1Smn=median(nonzeros(SRNcG1S)); NdifSmn=median(nonzeros(SRNcS)); 
NdifG1std=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifG1Sstd=std(SRNcG1S(SRNcG1S~=0))/sqrt(nnz(SRNcG1S)); NdifSstd=std(SRNcS(SRNcS~=0))/sqrt(nnz(SRNcS)); 
CdifG1mn=median(nonzeros(SRCyG1)); CdifG1Smn=median(nonzeros(SRCyG1S)); CdifSmn=median(nonzeros(SRCyS)); 
CdifG1std=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifG1Sstd=std(SRCyG1S(SRCyG1S~=0))/sqrt(nnz(SRCyG1S)); CdifSstd=std(SRCyS(SRCyS~=0))/sqrt(nnz(SRCyS)); 
WdifG1mn=median(nonzeros(SRWG1)); WdifG1Smn=median(nonzeros(SRWG1S)); WdifSmn=median(nonzeros(SRWS)); 
WdifG1std=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifG1Sstd=std(SRWS(SRWG1S~=0))/sqrt(nnz(SRWG1S)); WdifSstd=std(SRWS(SRWS~=0))/sqrt(nnz(SRWS)); 

Ny=[WdifG1mn,WdifG1Smn,WdifSmn;NdifG1mn,NdifG1Smn,NdifSmn;CdifG1mn,CdifG1Smn,CdifSmn];
% Ny=[NdifG1mn,NdifSmn,NdifG2mn;CdifG1mn,CdifSmn,CdifG2mn;WdifG1mn,WdifSmn,WdifG2mn];
Nerr=[WdifG1std,WdifG1Sstd,WdifSstd;NdifG1std,NdifG1Sstd,NdifSstd;CdifG1std,CdifG1Sstd,CdifSstd];
bar(Ny)
set(gca, 'XTickLabel', {'Whole cell','Nucleus','Cytoplasm'});
hold on
ngroups = size(Ny, 1);
nbars = size(Ny, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
ydata{1,1}=nonzeros(SRWG1); ydata{2,1}=nonzeros(SRNcG1); ydata{3,1}=nonzeros(SRCyG1);
ydata{1,2}=nonzeros(SRWG1S); ydata{2,2}=nonzeros(SRNcG1S); ydata{3,2}=nonzeros(SRCyG1S);
ydata{1,3}=nonzeros(SRWS); ydata{2,3}=nonzeros(SRNcS); ydata{3,3}=nonzeros(SRCyS);


for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Ny(:,i), Nerr(:,i), '.');
    for we=1:3
        [r, c] = size(ydata{we,i}'); xdata = repmat(1:c, r, 1);
        scatter(x(we).* ones(length(ydata{we,i}'),1), ydata{we,i}', 'k.', 'jitter','on', 'jitterAmount', 0.05);
    end
end
    
hold off
set(gcf,'color','w');
ylabel('Specific growth rate (/hr)');

set(gcf,'Color','w'); ylim([-0.4 0.4]);
ylabel('Specific growth rate (/hr)');

% S7.b: cellcycle phase growth over time in cell, nucleus, cytoplasm
load('F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\MassGenResults_rev136\SegmentedWSs\NCWsgr_rev1.mat')
plot(sgrN0(5:21,1)); hold on;
plot(sgrC0(5:21,1));
plot(sgrW0(5:21,1))

%S&.c: RFP-GFP FL intensity over cell cycle phases
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
fdirM='F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\MassGenResults_rev136\SegmentedWSs\';
fdir ='F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\Trial1\';
cellno=3;
for dd=6:22
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM],'BWNucl');
%     label(BWNucl);
    Lnls = BWNucl(QPIa:QPIb,QPIc:QPId);
    fnameg=sprintf('GFP120X_%d_frame_%d.tif',cellno,330+((dd-6)*60));
    GFP=imread([fdir fnameg]); GFP=GFP(:,1:FLd);
    GFP2 = single(imresize(GFP,[278,360]));
    GFP_bg = (sum(sum((1-Lnls).*GFP2)))/(sum(sum(1-Lnls)));
%     figure(1); imagesc(GFP); pause();
    IntG(dd-5)=sum(sum((GFP2./GFP_bg).*Lnls))/sum(sum(Lnls));
%     GFP3 = GFP2.*Lnls;
%     IntG(dd-5)=max(GFP3(:))/GFP_bg;
    fnamer=sprintf('RFP120X_%d_frame_%d.tif',cellno,330+((dd-6)*60));
    RFP=imread([fdir fnamer]); RFP=RFP(:,1:FLd); 
    RFP2 = single(imresize(RFP,[278,360]));
%     figure(2); imagesc(RFP); pause();
    RFP_bg = (sum(sum((1-Lnls).*RFP2)))/(sum(sum(1-Lnls)));
%     RFP3= RFP2.*Lnls;
%     IntR(dd-5)=max(RFP3(:))/RFP_bg;
%     figure(1); imagesc(GFP); pause();
    IntR(dd-5)=sum(sum((RFP2./RFP_bg).*Lnls))/sum(sum(Lnls));
end
plot(IntG/IntG(1)); hold on; plot(IntR/IntR(17));

%% figure S8: copy from figure 2, 3 and S5 with nuclear boundary
% figure S8c : S phase RPE

clear all; clc;
fdir='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdir1='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
mm=2; ii=14;

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
clf; fname=sprintf('WS%d_cell%d.mat',ii,mm);
load([fdir fname],'GC','Abkg_mass','Time','fstart','Abkg_stored2','MassAm','sz','xcg','kk'); 
for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
    for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
        xData=Time(1:kk+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
        yData=reshape(Abkg_mass(pp,qq,1:kk+1),1,kk+1);
        Indt=find(yData(:)==0,1);
        if isempty(Indt)
            Indt=kk+1;
            fitData1 = polyfit(xData(1:Indt), medfilt1(yData(1:Indt),5), 1);
            Massavg(pp,qq)=mean(yData(1:Indt));
            GC(pp,qq)=fitData1(1);
        else
            GC(pp,qq)=0;
        end
    end
end
SGf=zeros(512,512);
SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.001;
SGf = SGf.*SM;
DD=Abkg_stored2(:,:,1);
c = kmeans(DD(:), 3, 'MaxIter', 10000);
CC=reshape(c,512,512);
BWs1=CC~=mode(CC(:));
BWdfill = imfill(BWs1,8, 'holes');
seD = strel('diamond',1);
%BWfinal = imerode(BWnobord,seD);
BWfinal = imdilate(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWfinal=imerode(BWfinal,[se0 se90]);
 BWfinal = bwareaopen(BWfinal, 5000);
 se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);
fnameg=sprintf('mAG120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm);
%             fnameg=sprintf('GFP120X_%d_frame_%d.tif',mm,fstart-mm+1);
GFP=imread([fdir1 fnameg]); 
fnamer=sprintf('mKO2120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm); 
%             fnamer=sprintf('DAPI120X_%d_frame_%d.tif',mm,fstart-mm+1);
RFP=imread([fdir1 fnamer]);
RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
image=zeros(FLb,FLd,3);
mycolors1=zeros(2001,3);
mycolors1(1:2001,2)=0:0.0005:1;
%     maxg=single(min(GFP(:))+1000);
ming=single(min(GFP(:))+850);
maxg=single(max(GFP(:)))+150; 
%     maxg=2000; 
minr=single(min(RFP(:)))+400;
maxr=single(max(RFP(:)))+150; 
GFPmask1=GFP>ming;
GFPmask2=GFP<maxg;
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
image(:,:,2)=GFP2;

RFPmask1=RFP>minr;
RFPmask2=RFP<maxr;
RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,1)=RFP2;
DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
Image=rgb2gray(image);
BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
sz=size(BWn);
c = kmeans(BWn(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
jk=mode(CC(:));
BWsn=CC~=jk;
%     BWsn = imbinarize(BWn,0.15);
seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
BWfinaln = imerode(BWsn,seD);
BWfinaln = bwareaopen(BWfinaln, 200);
Ln1 = imclearborder(BWfinaln, 4);
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ln1), kernel, 'same');
binaryImage = blurryImage > 0.5; 
se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
%     binaryImage = imdilate(binaryImage, [se90 se0]);
%     Lnls=imresize(binaryImage,[356,458]);
Lnls=imresize(BWfinaln,[QPIb-QPIa+1,QPId-QPIc+1]);
BWNucl=zeros(512,512);
BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls; % before FL magnifier
%             BWNucl(131:486,55:512)=Lnls; % after FL magnifier
BWdilt=imdilate(BWNucl,[se0 se90]);
BWdilt= imerode(BWdilt,[se0 se90]);
figure(1);
imoverlay(Abkg_mass(:,:,1),SGf.*60,[-0.06, 0.06],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWdilt,'noholes');
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)


%% figure S8d: G2 phase

clear all; clc;
fdir='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdir1='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
mm=2; ii=18;

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
clf; fname=sprintf('WS%d_cell%d.mat',ii,mm);
load([fdir fname],'GC','Abkg_mass','Time','fstart','Abkg_stored2','MassAm','sz','xcg','kk'); 
for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
    for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
        xData=Time(1:kk+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
        yData=reshape(Abkg_mass(pp,qq,1:kk+1),1,kk+1);
        Indt=find(yData(:)==0,1);
        if isempty(Indt)
            Indt=kk+1;
            fitData1 = polyfit(xData(1:Indt), medfilt1(yData(1:Indt),5), 1);
            Massavg(pp,qq)=mean(yData(1:Indt));
            GC(pp,qq)=fitData1(1);
        else
            GC(pp,qq)=0;
        end
    end
end
SGf=zeros(512,512);
SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.001;
SGf = SGf.*SM;
DD=Abkg_stored2(:,:,1);
c = kmeans(DD(:), 3, 'MaxIter', 10000);
CC=reshape(c,512,512);
BWs1=CC~=mode(CC(:));
BWdfill = imfill(BWs1,8, 'holes');
seD = strel('diamond',1);
%BWfinal = imerode(BWnobord,seD);
BWfinal = imdilate(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWfinal=imerode(BWfinal,[se0 se90]);
 BWfinal = bwareaopen(BWfinal, 5000);
 se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);
fnameg=sprintf('mAG120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm);
%             fnameg=sprintf('GFP120X_%d_frame_%d.tif',mm,fstart-mm+1);
GFP=imread([fdir1 fnameg]); 
fnamer=sprintf('mKO2120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm); 
%             fnamer=sprintf('DAPI120X_%d_frame_%d.tif',mm,fstart-mm+1);
RFP=imread([fdir1 fnamer]);
RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
image=zeros(FLb,FLd,3);
mycolors1=zeros(2001,3);
mycolors1(1:2001,2)=0:0.0005:1;
%     maxg=single(min(GFP(:))+1000);
ming=single(min(GFP(:))+850);
maxg=single(max(GFP(:)))+150; 
%     maxg=2000; 
minr=single(min(RFP(:)))+400;
maxr=single(max(RFP(:)))+150; 
GFPmask1=GFP>ming;
GFPmask2=GFP<maxg;
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
image(:,:,2)=GFP2;

RFPmask1=RFP>minr;
RFPmask2=RFP<maxr;
RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,1)=RFP2;
DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
Image=rgb2gray(image);
BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
sz=size(BWn);
c = kmeans(BWn(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
jk=mode(CC(:));
BWsn=CC~=jk;
%     BWsn = imbinarize(BWn,0.15);
seD = strel('diamond',2);
%BWfinal = imerode(BWnobord,seD);
BWfinaln = imerode(BWsn,seD);
BWfinaln = bwareaopen(BWfinaln, 200);
Ln1 = imclearborder(BWfinaln, 4);
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ln1), kernel, 'same');
binaryImage = blurryImage > 0.5; 
se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
%     binaryImage = imdilate(binaryImage, [se90 se0]);
%     Lnls=imresize(binaryImage,[356,458]);
Lnls=imresize(BWfinaln,[QPIb-QPIa+1,QPId-QPIc+1]);
BWNucl=zeros(512,512);
BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls; % before FL magnifier
%             BWNucl(131:486,55:512)=Lnls; % after FL magnifier
BWdilt=imdilate(BWNucl,[se0 se90]);
BWdilt= imerode(BWdilt,[se0 se90]);
figure(1);
imoverlay(Abkg_mass(:,:,1),SGf.*60,[-0.06, 0.06],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWdilt,'noholes');
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)


%% figure S8g: S phase MCF7

% load('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\MassGenTrackingResults_rev136_Trial1\cell7\MassGenResults_rev13_6_trial2\WS20_cell7.mat')
% [BWfinal,SS] = imagebackground_poly4_bead(D_stored2(:,:,1));
% SGf=zeros(512,512); SGf = imfilter((GC.*60), fspecial('gaussian', [50 50], 1));
% SM = SS(5:508,5:508)>-0.05; SGf = SGf.*SM; SGf(SGf==0) = NaN;
% figure(5);
% imoverlay(SS(5:508,5:508),SGf,[-0.03, 0.03],[-0.05 0.35],parula, 0.2, gca); 
% colormap(map); ylim([190 490]); xlim([104 404]);
% GFP=imread('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\Trial1\GFP120X_7_frame_1170.tif');
% RFP=imread('F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\Trial1\RFP120X_7_frame_1170.tif');

clear all; clc;
fdir='F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\MassGenResults_rev136\SegmentedWSs\';
mm=4; ii=14;
clf; fname=sprintf('WS%d_cell%d.mat',ii,mm);
load([fdir fname],'GC','Abkg_mass','Time','fstart','BWNucl','sz','xcg','kk'); 
for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
    for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
        xData=Time(1:kk+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
        yData=reshape(Abkg_mass(pp,qq,1:kk+1),1,kk+1);
        Indt=find(yData(:)==0,1);
        if isempty(Indt)
            Indt=kk+1;
            fitData1 = polyfit(xData(1:Indt), medfilt1(yData(1:Indt),5), 1);
            Massavg(pp,qq)=mean(yData(1:Indt));
            GC(pp,qq)=fitData1(1);
        else
            GC(pp,qq)=0;
        end
    end
end
SGf=zeros(512,512);
SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.001;
SGf = SGf.*SM;

figure(1);
imoverlay(Abkg_mass(:,:,1),SGf.*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWNucl,'noholes');
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)

%% figure S12

% Doxorubicin drug impact on RPE
filenow{1}='S:\Data\Soorya\RPEDrugTest_2020\Pn10_Water2.3uL_18June2020\MassGenresults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{2}='K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{3}='K:\Data\Soorya\RPEDrugTest_2020\Pn4_Doxorubicin0.01um_30April2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{4}='S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin0.01uM_12June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{5}='K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{6}='S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin0.2uM_13June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{7}='K:\Data\Soorya\RPEDrugTest_2020\Pn5_Doxorubicin4um_1May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{8}='S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin4uM_13June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';

for filew=1:4
    load([filenow{(filew-1)+1}]);
    for yy=1:numLoc
        mass=MassWhole(yy,:); Time=TimeWhole(yy,:);
        if ErrPoint(yy)>120
            fitData=polyfit(Time(1:ErrPoint(yy)),mass(1:ErrPoint(yy)),1);
            sgrfilew(filew,yy)=(fitData(1)/mean(mass(1:ErrPoint(yy))))*60;
        else
            sgrfilew(filew,yy)=0;
        end
    end
    load([filenow{(filew*2)}]);
    for zz=yy+1:yy+numLoc
        mass=MassWhole(zz-yy,:); Time=TimeWhole(zz-yy,:);
        if ErrPoint(zz-yy)>120
            fitData=polyfit(Time(1:ErrPoint(zz-yy)),mass(1:ErrPoint(zz-yy)),1);
            sgrfilew(filew,zz)=(fitData(1)/mean(mass(1:ErrPoint(zz-yy))))*60;
        else
            sgrfilew(filew,zz)=0;
        end
    end
end
sgrfilew(isnan(sgrfilew))=0;
%figure(2);
for xx=1:4
    sgrmN(xx)=mean(nonzeros(sgrfilew(xx,:)));
    sgrSt(xx)=std(nonzeros(sgrfilew(xx,:)));
    sgrVl{xx}=nonzeros(sgrfilew(xx,:));
end
% bar(sgrmN);
% set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2','Dox 4'});
% hold on;
% errorbar([1,2,3,4],sgrmN,sgrSt,sgrSt,'.','Color','k');
% set(gcf,'Color','w');
figure(1); violin(sgrVl);
set(gcf,'color','w'); ylim([-0.12 0.12]);

% Homoharringtonine drug impact on RPE
filenow{1}='K:\Data\Soorya\RPEDrugTest_2020\Pn6_DMSO2.1uL_4May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{2}='S:\Data\Soorya\RPEDrugTest_2020\Pn10_DMSO2.1uL_18June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{3}='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.036um_3May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{4}='T:\Data\Soorya\RPEDrugTest2020\Pn7_Homoharringtonine0.036uM_13Jun2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{5}='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.36um_3May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{6}='S:\Data\Soorya\RPEDrugTest_2020\Pn8_Homoharringtonine0.36uM_14Jun2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{7}='K:\Data\Soorya\RPEDrugTest_2020\Pn3Homoharringtonine_21March2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{8}='S:\Data\Soorya\RPEDrugTest_2020\Pn8_Homoharringtonine3.6uM_14Jun2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';

for filew=1:4
    load([filenow{(filew-1)+1}]);
    for yy=1:numLoc
        mass=MassWhole(yy,:); Time=TimeWhole(yy,:);
        if ErrPoint(yy)>120
            fitData=polyfit(Time(1:ErrPoint(yy)),mass(1:ErrPoint(yy)),1);
            sgrfilew(filew,yy)=(fitData(1)/mean(mass(1:ErrPoint(yy))))*60;
        else
            sgrfilew(filew,yy)=0;
        end
    end
    load([filenow{(filew*2)}]);
    for zz=yy+1:yy+numLoc
        mass=MassWhole(zz-yy,:); Time=TimeWhole(zz-yy,:);
        if ErrPoint(zz-yy)>120
            fitData=polyfit(Time(1:ErrPoint(zz-yy)),mass(1:ErrPoint(zz-yy)),1);
            sgrfilew(filew,zz)=(fitData(1)/mean(mass(1:ErrPoint(zz-yy))))*60;
        else
            sgrfilew(filew,zz)=0;
        end
    end
end
sgrfilew(isnan(sgrfilew))=0;
% figure(2);
for xx=1:4
    sgrmN(xx)=mean(nonzeros(sgrfilew(xx,:)));
    sgrSt(xx)=std(nonzeros(sgrfilew(xx,:)));
    sgrVl{xx}=nonzeros(sgrfilew(xx,:));
end
% bar(sgrmN);
% set(gca, 'XTickLabel', {'DMSO','Homo 0.036','Homo 0.36','Homo 3.6'});
% hold on;
% errorbar([1,2,3,4],sgrmN,sgrSt,sgrSt,'.','Color','k');
% set(gcf,'Color','w');
figure(2); violin(sgrVl);
set(gcf,'color','w'); ylim([-0.12 0.12]);

% Cycloheximide drug impact on RPE
filenow{1}='K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{2}='S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{3}='K:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.1um_21May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{4}='T:\Data\Soorya\RPEDrugTest2020\Pn7_Cycloheximide0.1uM_11June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{5}='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide0.53um_20May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{6}='S:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.53uM_11June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{7}='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide2.88um_20May2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';
filenow{8}='T:\Data\Soorya\RPEDrugTest2020\Pn7_Cycloheximide2.88uM_11June2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat';

for filew=1:4
    load([filenow{(filew-1)+1}]);
    for yy=1:numLoc
        mass=MassWhole(yy,:); Time=TimeWhole(yy,:);
        if ErrPoint(yy)>120
            fitData=polyfit(Time(1:ErrPoint(yy)),mass(1:ErrPoint(yy)),1);
            sgrfilew(filew,yy)=(fitData(1)/mean(mass(1:ErrPoint(yy))))*60;
        else
            sgrfilew(filew,yy)=0;
        end
    end
    load([filenow{(filew*2)}]);
    for zz=yy+1:yy+numLoc
        mass=MassWhole(zz-yy,:); Time=TimeWhole(zz-yy,:);
        if ErrPoint(zz-yy)>120
            fitData=polyfit(Time(1:ErrPoint(zz-yy)),mass(1:ErrPoint(zz-yy)),1);
            sgrfilew(filew,zz)=(fitData(1)/mean(mass(1:ErrPoint(zz-yy))))*60;
        else
            sgrfilew(filew,zz)=0;
        end
    end
end
sgrfilew(isnan(sgrfilew))=0;
% figure(2);
for xx=1:4
    sgrmN(xx)=mean(nonzeros(sgrfilew(xx,:)));
    sgrSt(xx)=std(nonzeros(sgrfilew(xx,:)));
    sgrVl{xx}=nonzeros(sgrfilew(xx,:));
end
% bar(sgrmN);
% set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53','Cyclo 2.88'});
% hold on;
% errorbar([1,2,3,4],sgrmN,sgrSt,sgrSt,'.','Color','k');
% set(gcf,'Color','w');
figure(3); violin(sgrVl);
set(gcf,'color','w'); ylim([-0.12 0.12]);


% figure S12d :  impact of dox on RPE over time

% temporal effect of doxorubicin, Homo and Cyclo
load('T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Mar2022\WS_S12de.mat'); % Dox
% load('T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Mar2022\WS_S12fg.mat'); % Homo
% load('T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Mar2022\WS_S12hi.mat'); % Cyclo

hr1=sgr1(:,1);
hr2=sgr1(:,2);
hr3=sgr1(:,3);
hr4=sgr1(:,4);
figure(1);
bar([median(hr1(hr1~=0)),median(hr2(hr2~=0)),median(hr3(hr3~=0)),median(hr4(hr4~=0))]);
set(gca, 'XTickLabel', {'First','Second','Third','Fourth'});
hold on;
errorbar([1,2,3,4],[median(hr1(hr1~=0)),median(hr2(hr2~=0)),median(hr3(hr3~=0)),median(hr4(hr4~=0))],...
    [std(hr1(hr1~=0))/sqrt(nnz(hr1)),std(hr2(hr2~=0))/sqrt(nnz(hr2)),std(hr3(hr3~=0))/sqrt(nnz(hr3)),std(hr4(hr4~=0))/sqrt(nnz(hr4))],...
    [std(hr1(hr1~=0))/sqrt(nnz(hr1)),std(hr2(hr2~=0))/sqrt(nnz(hr2)),std(hr3(hr3~=0))/sqrt(nnz(hr3)),std(hr4(hr4~=0))/sqrt(nnz(hr4))],...
    '.','Color','b');
set(gcf,'Color','w'); ylim([-0.1 0.1]);
ylabel('Specific growth rate (/hr)');

hr1=sgr2(:,1);
hr2=sgr2(:,2);
hr3=sgr2(:,3);
hr4=sgr2(:,4);
figure(2);
bar([median(hr1(hr1~=0)),median(hr2(hr2~=0)),median(hr3(hr3~=0)),median(hr4(hr4~=0))]);
set(gca, 'XTickLabel', {'First','Second','Third','Fourth'});
hold on;
errorbar([1,2,3,4],[median(hr1(hr1~=0)),median(hr2(hr2~=0)),median(hr3(hr3~=0)),median(hr4(hr4~=0))],...
    [std(hr1(hr1~=0))/sqrt(nnz(hr1)),std(hr2(hr2~=0))/sqrt(nnz(hr2)),std(hr3(hr3~=0))/sqrt(nnz(hr3)),std(hr4(hr4~=0))/sqrt(nnz(hr4))],...
    [std(hr1(hr1~=0))/sqrt(nnz(hr1)),std(hr2(hr2~=0))/sqrt(nnz(hr2)),std(hr3(hr3~=0))/sqrt(nnz(hr3)),std(hr4(hr4~=0))/sqrt(nnz(hr4))],...
    '.','Color','b');
set(gcf,'Color','w'); ylim([-0.1 0.1]);
ylabel('Specific growth rate (/hr)');

%% figure S13 : Doxorubicin drug impact on RPE cells

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

% figure S13 a, b, c, d
fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn4_Doxorubicin0.01um_30April2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn4_Doxorubicin0.01um_30April2020\Trial1\';
mm=29; % 23
ii=4;

% for Dox 0.2
fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\Trial1\';
mm=3; % 23
ii=2;

% for Dox 4
fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn5_Doxorubicin4um_1May2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn5_Doxorubicin4um_1May2020\Trial1\';
mm=5; % 23
ii=1;

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
clf; fname=sprintf('WS%d_cell%d.mat',ii,mm);
if exist([fdir fname])~=0
    load([fdir fname],'GC','Abkg_mass','Time','fstart','Abkg_stored2','MassAm','sz','xcg','kk'); 
    for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
        for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
            xData=Time(1:kk+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
            yData=reshape(Abkg_mass(pp,qq,1:kk+1),1,kk+1);
            Indt=find(yData(:)==0,1);
            if isempty(Indt)
                Indt=kk+1;
                fitData1 = polyfit(xData(1:Indt), medfilt1(yData(1:Indt),5), 1);
                Massavg(pp,qq)=mean(yData(1:Indt));
                GC(pp,qq)=fitData1(1);
            else
                GC(pp,qq)=0;
            end
        end
    end
    SGf=zeros(512,512);
    SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    SM = Abkg_mass(:,:,1)>-0.001;
    SGf = SGf.*SM;
    DD=Abkg_stored2(:,:,1);
    c = kmeans(DD(:), 3, 'MaxIter', 10000);
    CC=reshape(c,512,512);
    BWs1=CC~=mode(CC(:));
    BWdfill = imfill(BWs1,8, 'holes');
    seD = strel('diamond',1);
    %BWfinal = imerode(BWnobord,seD);
    BWfinal = imdilate(BWdfill,seD);
    BWfinal = bwareaopen(BWfinal, 2000);
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    BWfinal=imerode(BWfinal,[se0 se90]);
     BWfinal = bwareaopen(BWfinal, 5000);
     se90 = strel('line', 8, 90);
    se0 = strel('line', 8, 0);
    BWfinal=imdilate(BWfinal,[se0 se90]);
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm);
%             fnameg=sprintf('GFP120X_%d_frame_%d.tif',mm,fstart-mm+1);
    GFP=imread([fdir1 fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm); 
%             fnamer=sprintf('DAPI120X_%d_frame_%d.tif',mm,fstart-mm+1);
    RFP=imread([fdir1 fnamer]);
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    image=zeros(FLb,FLd,3);
    mycolors1=zeros(2001,3);
    mycolors1(1:2001,2)=0:0.0005:1;
    %     maxg=single(min(GFP(:))+1000);
    ming=single(min(GFP(:))+850);
    maxg=single(max(GFP(:)))+150; 
    %     maxg=2000; 
    minr=single(min(RFP(:)))+400;
    maxr=single(max(RFP(:)))-250; 
    GFPmask1=GFP>ming;
    GFPmask2=GFP<maxg;
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
    GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
    image(:,:,2)=GFP2;

    RFPmask1=RFP>minr;
    RFPmask2=RFP<maxr;
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
    RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
    mycolors2=zeros(2001,3);
    mycolors2(1:2001,1)=0:0.0005:1;
    image(:,:,1)=RFP2;
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    Image=rgb2gray(image);
    BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
    sz=size(BWn);
    c = kmeans(BWn(:), 2, 'MaxIter', 10000);
    CC=reshape(c,sz(1),sz(2));
    jk=mode(CC(:));
    BWsn=CC~=jk;
    %     BWsn = imbinarize(BWn,0.15);
    seD = strel('diamond',2);
    %BWfinal = imerode(BWnobord,seD);
    BWfinaln = imerode(BWsn,seD);
    BWfinaln = bwareaopen(BWfinaln, 200);
    Ln1 = imclearborder(BWfinaln, 4);
    windowSize = 51;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(Ln1), kernel, 'same');
    binaryImage = blurryImage > 0.5; 
    se90 = strel('line', 30, 90);
    se0 = strel('line', 30, 0);
    %     binaryImage = imdilate(binaryImage, [se90 se0]);
    %     Lnls=imresize(binaryImage,[356,458]);
    Lnls=imresize(BWfinaln,[QPIb-QPIa+1,QPId-QPIc+1]);
    BWNucl=zeros(512,512);
    BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls; % before FL magnifier
%             BWNucl(131:486,55:512)=Lnls; % after FL magnifier
    BWdilt=imdilate(BWNucl,[se0 se90]);
    BWCyto=BWfinal.*(1-BWNucl).*(BWdilt);
    figure(2);
    imoverlay(Abkg_mass(135:430,117:416,1),SGf(135:430,117:416).*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
    colormap(temp); hold on;
    [Bn,Lnlls] = bwboundaries(BWNucl(135:430,117:416),'noholes');
   for ee=1:length(Bn)
        boundary = Bn{ee};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
    end
%             DImage2=imresize(DImage,[356,458,3]);
    %             [B,Llls] = bwboundaries(BWfinal(135:430,117:416),'noholes');
%             hold on
%             for ee=1:length(B)
%                 boundary = B{ee};
%                 plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
%             end
%             DImage2=imresize(DImage,[356,458,3]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    figure(1)
    imshow(DImage3(135:430,117:416,:)); hold on;
%             for ee=1:length(Bn)
%                 boundary = Bn{ee};
%                 plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
%             end
%             [B,Llls] = bwboundaries(BWfinal(135:430,117:416),'noholes');
%             hold on
%             for ee=1:length(B)
%                 boundary = B{ee};
%                 plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
%             endDImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=220; ybase=240; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textdown=sprintf('10 µm');
    text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
    pause(2);
end

% figure S13 e: overall growth rate at different Dox dose
load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1=sgrN0.*InDepth; sgrC0_1=sgrC0.*InDepth; sgrW0_1=sgrW0.*InDepth; CellPhase_1=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn10_Water2.3uL_18June2020\MassGenresults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1(31:60,:)=sgrN0.*InDepth; sgrC0_1(31:60,:)=sgrC0.*InDepth; sgrW0_1(31:60,:)=sgrW0.*InDepth; CellPhase_1(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn4_Doxorubicin0.01um_30April2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2=sgrN0.*InDepth; sgrC0_2=sgrC0.*InDepth; sgrW0_2=sgrW0.*InDepth; CellPhase_2=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin0.01uM_12June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2(31:60,:)=sgrN0.*InDepth; sgrC0_2(31:60,:)=sgrC0.*InDepth; sgrW0_2(31:60,:)=sgrW0.*InDepth; CellPhase_2(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3=sgrN0.*InDepth; sgrC0_3=sgrC0.*InDepth; sgrW0_3=sgrW0.*InDepth; CellPhase_3=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin0.2uM_13June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3(31:60,:)=sgrN0.*InDepth; sgrC0_3(31:60,:)=sgrC0.*InDepth; sgrW0_3(31:60,:)=sgrW0.*InDepth; CellPhase_3(31:60,:)=CellPhase.*InDepth;
% load('K:\Data\Soorya\RPEDrugTest_2020\Pn5_Doxorubicin4um_1May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev4.mat');
% sgrN0_4=sgrN0; sgrC0_4=sgrC0; sgrW0_4=sgrW0; CellPhase_4=CellPhase;
% load('S:\Data\Soorya\RPEDrugTest_2020\Pn7_Doxorubicin4uM_13June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev4.mat');
% sgrN0_4(31:60,:)=sgrN0; sgrC0_4(31:60,:)=sgrC0; sgrW0_4(31:60,:)=sgrW0; CellPhase_4(31:60,:)=CellPhase;

group1=[mean(sgrN0_1(sgrN0_1~=0)),mean(sgrN0_2(sgrN0_2~=0)),mean(sgrN0_3(sgrN0_3~=0))];
group2=[mean(sgrC0_1(sgrC0_1~=0)),mean(sgrC0_2(sgrC0_2~=0)),mean(sgrC0_3(sgrC0_3~=0))];
group3=[mean(sgrW0_1(sgrW0_1~=0)),mean(sgrW0_2(sgrW0_2~=0)),mean(sgrW0_3(sgrW0_3~=0))];
Egroup1=[std(sgrN0_1(sgrN0_1~=0))/sqrt(nnz(sgrN0_1)),std(sgrN0_2(sgrN0_2~=0))/sqrt(nnz(sgrN0_2)),std(sgrN0_3(sgrN0_3~=0))/sqrt(nnz(sgrN0_3))];
Egroup2=[std(sgrC0_1(sgrC0_1~=0))/sqrt(nnz(sgrC0_1)),std(sgrC0_2(sgrC0_2~=0))/sqrt(nnz(sgrC0_2)),std(sgrC0_3(sgrC0_3~=0))/sqrt(nnz(sgrC0_3))];
Egroup3=[std(sgrW0_1(sgrW0_1~=0))/sqrt(nnz(sgrW0_1)),std(sgrW0_2(sgrW0_2~=0))/sqrt(nnz(sgrW0_2)),std(sgrW0_3(sgrW0_3~=0))/sqrt(nnz(sgrW0_3))];
y=[group1;group2;group3];
err=[Egroup1;Egroup2;Egroup3];
for ww=1:3
    errorbar(y(ww,:),err(ww,:));
    hold on;
end
set(gcf,'color','w'); xlim([0 4]);
ylabel('Specific growth rate (/hr)'); xticks([1 2 3])
set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2'});
ylim([-0.15 0.15]);


% figure S13 f,g,h: Nuclear sgr in G1, S, G2 phases 

G1WS=(CellPhase_1==1); G1SWS=(CellPhase_1==2); SWSG2=(CellPhase_1==3 | CellPhase_1==4);
SRNcG1=G1WS.*sgrN0_1; SRNcG1S=G1SWS.*sgrN0_1; SRNcSG2=SWSG2.*sgrN0_1;
NdifG1mn_1=median(nonzeros(SRNcG1)); NdifG1Smn_1=median(nonzeros(SRNcG1S)); NdifSG2mn_1=median(nonzeros(SRNcSG2)); 
NdifG1std_1=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifG1Sstd_1=std(SRNcG1S(SRNcG1S~=0))/sqrt(nnz(SRNcG1S)); NdifSG2std_1=std(SRNcSG2(SRNcSG2~=0))/sqrt(nnz(SRNcSG2)); 

SRCyG1=G1WS.*sgrC0_1; SRCyG1S=G1SWS.*sgrC0_1; SRCySG2=SWSG2.*sgrC0_1; 
CdifG1mn_1=median(nonzeros(SRCyG1)); CdifG1Smn_1=median(nonzeros(SRCyG1S)); CdifSG2mn_1=median(nonzeros(SRCySG2)); 
CdifG1std_1=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifG1Sstd_1=std(SRCyG1S(SRCyG1S~=0))/sqrt(nnz(SRCyG1S)); CdifSG2std_1=std(SRCySG2(SRCySG2~=0))/sqrt(nnz(SRCySG2)); 

SRWG1=G1WS.*sgrW0_1; SRWG1S=G1SWS.*sgrW0_1; SRWSG2=SWSG2.*sgrW0_1;
WdifG1mn_1=median(nonzeros(SRWG1)); WdifG1Smn_1=median(nonzeros(SRWG1S)); WdifSG2mn_1=median(nonzeros(SRWSG2)); 
WdifG1std_1=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifG1Sstd_1=std(SRWG1S(SRWG1S~=0))/sqrt(nnz(SRWG1S)); WdifSG2std_1=std(SRWSG2(SRWSG2~=0))/sqrt(nnz(SRWSG2)); 

G1WS=(CellPhase_2==1); G1SWS=(CellPhase_2==2); SWSG2=(CellPhase_2==3 | CellPhase_2==4); 

SRNcG1=G1WS.*sgrN0_2; SRNcG1S=G1SWS.*sgrN0_2; SRNcSG2=SWSG2.*sgrN0_2; 
NdifG1mn_2=median(nonzeros(SRNcG1)); NdifG1Smn_2=median(nonzeros(SRNcG1S)); NdifSG2mn_2=median(nonzeros(SRNcSG2)); 
NdifG1std_2=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifG1Sstd_2=std(SRNcG1S(SRNcG1S~=0))/sqrt(nnz(SRNcG1S)); NdifSG2std_2=std(SRNcSG2(SRNcSG2~=0))/sqrt(nnz(SRNcSG2)); 


SRCyG1=G1WS.*sgrC0_2; SRCyG1S=G1SWS.*sgrC0_2; SRCySG2=SWSG2.*sgrC0_2; 
CdifG1mn_2=median(nonzeros(SRCyG1)); CdifG1Smn_2=median(nonzeros(SRCyG1S)); CdifSG2mn_2=median(nonzeros(SRCySG2)); 
CdifG1std_2=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifG1Sstd_2=std(SRCyG1S(SRCyG1S~=0))/sqrt(nnz(SRCyG1S)); CdifSG2std_2=std(SRCySG2(SRCySG2~=0))/sqrt(nnz(SRCySG2)); 

SRWG1=G1WS.*sgrW0_2; SRWG1S=G1SWS.*sgrW0_2; SRWSG2=SWSG2.*sgrW0_2; 
WdifG1mn_2=median(nonzeros(SRWG1)); WdifG1Smn_2=median(nonzeros(SRWG1S)); WdifSG2mn_2=median(nonzeros(SRWSG2));
WdifG1std_2=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifG1Sstd_2=std(SRWG1S(SRWG1S~=0))/sqrt(nnz(SRWG1S)); WdifSG2std_2=std(SRWSG2(SRWSG2~=0))/sqrt(nnz(SRWSG2)); 

G1WS=(CellPhase_3==1); G1SWS=(CellPhase_3==2); SWSG2=(CellPhase_3==3 | CellPhase_3==4); 

SRNcG1=G1WS.*sgrN0_3; SRNcG1S=G1SWS.*sgrN0_3; SRNcSG2=SWSG2.*sgrN0_3; 
NdifG1mn_3=median(nonzeros(SRNcG1)); NdifG1Smn_3=median(nonzeros(SRNcG1S)); NdifSG2mn_3=median(nonzeros(SRNcSG2));
NdifG1std_3=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifG1Sstd_3=std(SRNcG1S(SRNcG1S~=0))/sqrt(nnz(SRNcG1S)); NdifSG2std_3=std(SRNcSG2(SRNcSG2~=0))/sqrt(nnz(SRNcSG2)); 

SRCyG1=G1WS.*sgrC0_3; SRCyG1S=G1SWS.*sgrC0_3; SRCySG2=SWSG2.*sgrC0_3; 
CdifG1mn_3=median(nonzeros(SRCyG1)); CdifG1Smn_3=median(nonzeros(SRCyG1S)); CdifSG2mn_3=median(nonzeros(SRCySG2));
CdifG1std_3=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifG1Sstd_3=std(SRCyG1S(SRCyG1S~=0))/sqrt(nnz(SRCyG1S)); CdifSG2std_3=std(SRCySG2(SRCySG2~=0))/sqrt(nnz(SRCySG2)); 

SRWG1=G1WS.*sgrW0_3; SRWG1S=G1SWS.*sgrW0_3; SRWSG2=SWSG2.*sgrW0_3; 
WdifG1mn_3=median(nonzeros(SRWG1)); WdifG1Smn_3=median(nonzeros(SRWG1S)); WdifSG2mn_3=median(nonzeros(SRWSG2)); 
WdifG1std_3=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifG1Sstd_3=std(SRWG1S(SRWG1S~=0))/sqrt(nnz(SRWG1S)); WdifSG2std_3=std(SRWSG2(SRWSG2~=0))/sqrt(nnz(SRWSG2));

% G1WS=(CellPhase_4==1); SWS=(CellPhase_4==2); G2WS=(CellPhase_4==3);
% 
% SRNcG1=G1WS.*sgrN0_4; SRNcS=SWS.*sgrN0_4; SRNcG2=G2WS.*sgrN0_4;
% NdifG1mn_4=median(nonzeros(SRNcG1)); NdifSmn_4=median(nonzeros(SRNcS)); NdifG2mn_4=median(nonzeros(SRNcG2));
% NdifG1std_4=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifSstd_4=std(SRNcS(SRNcS~=0))/sqrt(nnz(SRNcS)); NdifG2std_4=std(SRNcG2(SRNcG2~=0))/sqrt(nnz(SRNcG2));
% 
% SRCyG1=G1WS.*sgrC0_4; SRCyS=SWS.*sgrC0_4; SRCyG2=G2WS.*sgrC0_4;
% CdifG1mn_4=median(nonzeros(SRCyG1)); CdifSmn_4=median(nonzeros(SRCyS)); CdifG2mn_4=median(nonzeros(SRCyG2));
% CdifG1std_4=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifSstd_4=std(SRCyS(SRCyS~=0))/sqrt(nnz(SRCyS)); CdifG2std_4=std(SRCyG2(SRCyG2~=0))/sqrt(nnz(SRCyG2));
% 
% SRWG1=G1WS.*sgrW0_4; SRWS=SWS.*sgrW0_4; SRWG2=G2WS.*sgrW0_4;
% WdifG1mn_4=median(nonzeros(SRWG1)); WdifSmn_4=median(nonzeros(SRWS)); WdifG2mn_4=median(nonzeros(SRWG2));
% WdifG1std_4=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifSstd_4=std(SRWS(SRWS~=0))/sqrt(nnz(SRWS)); WdifG2std_4=std(SRWG2(SRWG2~=0))/sqrt(nnz(SRWG2));

group0=[mean(sgrN0_1(sgrN0_1~=0)),mean(sgrN0_2(sgrN0_2~=0)),mean(sgrN0_3(sgrN0_3~=0))]; %mean(sgrN0_4(sgrN0_4~=0))];
Egroup0=[std(sgrN0_1(sgrN0_1~=0))/sqrt(nnz(sgrN0_1)),std(sgrN0_2(sgrN0_2~=0))/sqrt(nnz(sgrN0_2)),std(sgrN0_3(sgrN0_3~=0))/sqrt(nnz(sgrN0_3))]; %,std(sgrN0_4(sgrN0_4~=0))/sqrt(nnz(sgrN0_4))];
group1=[mean(sgrC0_1(sgrC0_1~=0)),mean(sgrC0_2(sgrC0_2~=0)),mean(sgrC0_3(sgrC0_3~=0))];
Egroup1=[std(sgrC0_1(sgrC0_1~=0))/sqrt(nnz(sgrC0_1)),std(sgrC0_2(sgrC0_2~=0))/sqrt(nnz(sgrC0_2)),std(sgrC0_3(sgrC0_3~=0))/sqrt(nnz(sgrC0_3))];
group2=[mean(sgrW0_1(sgrW0_1~=0)),mean(sgrW0_2(sgrW0_2~=0)),mean(sgrW0_3(sgrW0_3~=0))];
Egroup2=[std(sgrW0_1(sgrW0_1~=0))/sqrt(nnz(sgrW0_1)),std(sgrW0_2(sgrW0_2~=0))/sqrt(nnz(sgrW0_2)),std(sgrW0_3(sgrW0_3~=0))/sqrt(nnz(sgrW0_3))];

figure(1); y=[group0;group1;group2];
err=[Egroup0;Egroup1;Egroup2];
for ww=1:3
    errorbar(y(ww,:),err(ww,:));
    hold on;
end
set(gcf,'color','w'); ylim([-0.15 0.15]); xlim([0 4]);
ylabel('Specific growth rate (/hr)'); xticks([1 2 3]); title('All phase');
set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2'});


group1=[NdifG1mn_1,NdifG1mn_2,NdifG1mn_3]; %,NdifG1mn_4]; 
Egroup1=[NdifG1std_1, NdifG1std_2,NdifG1std_3];  %,NdifG1std_4];
group2=[CdifG1mn_1,CdifG1mn_2,CdifG1mn_3]; 
Egroup2=[CdifG1std_1, CdifG1std_2,CdifG1std_3];
group3=[WdifG1mn_1,WdifG1mn_2,WdifG1mn_3]; 
Egroup3=[WdifG1std_1, WdifG1std_2,WdifG1std_3];

figure(2); y=[group1;group2;group3];
err=[Egroup1;Egroup2;Egroup3];
for ww=1:3
    errorbar(y(ww,:),err(ww,:));
    hold on;
end
set(gcf,'color','w'); ylim([-0.15 0.15]); xlim([0 4]);
ylabel('Specific growth rate (/hr)'); xticks([1 2 3]); title('G1 phase');
set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2'});

group1=[NdifG1Smn_1,NdifG1Smn_2,NdifG1Smn_3]; %,NdifSmn_4]; 
Egroup1=[NdifG1Sstd_1, NdifG1Sstd_2,NdifG1Sstd_3]; %,NdifSstd_4];
group2=[CdifG1Smn_1,CdifG1Smn_2,CdifG1Smn_3]; 
Egroup2=[CdifG1Sstd_1, CdifG1Sstd_2,CdifG1Sstd_3];
group3=[WdifG1Smn_1,WdifG1Smn_2,WdifG1Smn_3]; 
Egroup3=[WdifG1Sstd_1, WdifG1Sstd_2,WdifG1Sstd_3];

figure(3); y=[group1;group2;group3];
err=[Egroup1;Egroup2;Egroup3];
for ww=1:3
    errorbar(y(ww,:),err(ww,:));
    hold on;
end
set(gcf,'color','w'); ylim([-0.15 0.15]); xlim([0 4]);
ylabel('Specific growth rate (/hr)'); xticks([1 2 3]); title('G1S phase');
set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2'});

group1=[NdifSG2mn_1,NdifSG2mn_2,NdifSG2mn_3]; %,NdifSmn_4]; 
Egroup1=[NdifSG2std_1, NdifSG2std_2,NdifSG2std_3]; %,NdifSstd_4];
group2=[CdifSG2mn_1,CdifSG2mn_2,CdifSG2mn_3]; 
Egroup2=[CdifSG2std_1, CdifSG2std_2,CdifSG2std_3];
group3=[WdifSG2mn_1,WdifSG2mn_2,WdifSG2mn_3]; 
Egroup3=[WdifSG2std_1, WdifSG2std_2,WdifSG2std_3];

figure(4); y=[group1;group2;group3];
err=[Egroup1;Egroup2;Egroup3];
for ww=1:3
    errorbar(y(ww,:),err(ww,:));
    hold on;
end
set(gcf,'color','w'); ylim([-0.15 0.15]); xlim([0 4]);
ylabel('Specific growth rate (/hr)'); xticks([1 2 3]); title('S/G2 phase');
set(gca, 'XTickLabel', {'Water','Dox 0.01','Dox 0.2'});

%% figure S14 :  Homoharringtonine drug impact on RPE cells

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

% for Homo 0.036
fdir='T:\Data\Soorya\RPEDrugTest2020\Pn7_Homoharringtonine0.036uM_13Jun2020\MassGenResults_rev136\Results_2hr\';
fdir1='T:\Data\Soorya\RPEDrugTest2020\Pn7_Homoharringtonine0.036uM_13Jun2020\Trial1\';
mm=15; ii=2;

% for Homo 0.36
fdir='S:\Data\Soorya\RPEDrugTest_2020\Pn8_Homoharringtonine0.36uM_14Jun2020\MassGenResults_rev136\Results_2hr\';
fdir1='S:\Data\Soorya\RPEDrugTest_2020\Pn8_Homoharringtonine0.36uM_14Jun2020\Trial1\';
mm=23; ii=2;

% for Homo 3.6
% fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn3Homoharringtonine_21March2020\MassGenResults_rev136\Results_2hr\';
% fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn3Homoharringtonine_21March2020\Trial2\';
% mm=8; ii=1;

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
clf; fname=sprintf('WS%d_cell%d.mat',ii,mm);
if exist([fdir fname])~=0
    load([fdir fname],'GC','Abkg_mass','Time','fstart','Abkg_stored2','MassAm','sz','xcg','kk'); 
    for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
        for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
            xData=Time(1:kk+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
            yData=reshape(Abkg_mass(pp,qq,1:kk+1),1,kk+1);
            Indt=find(yData(:)==0,1);
            if isempty(Indt)
                Indt=kk+1;
                fitData1 = polyfit(xData(1:Indt), medfilt1(yData(1:Indt),5), 1);
                Massavg(pp,qq)=mean(yData(1:Indt));
                GC(pp,qq)=fitData1(1);
            else
                GC(pp,qq)=0;
            end
        end
    end
    SGf=zeros(512,512);
    SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    SM = Abkg_mass(:,:,1)>-0.001;
    SGf = SGf.*SM;
    DD=Abkg_stored2(:,:,1);
    c = kmeans(DD(:), 3, 'MaxIter', 10000);
    CC=reshape(c,512,512);
    BWs1=CC~=mode(CC(:));
    BWdfill = imfill(BWs1,8, 'holes');
    seD = strel('diamond',1);
    %BWfinal = imerode(BWnobord,seD);
    BWfinal = imdilate(BWdfill,seD);
    BWfinal = bwareaopen(BWfinal, 2000);
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    BWfinal=imerode(BWfinal,[se0 se90]);
     BWfinal = bwareaopen(BWfinal, 5000);
     se90 = strel('line', 8, 90);
    se0 = strel('line', 8, 0);
    BWfinal=imdilate(BWfinal,[se0 se90]);
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm);
%             fnameg=sprintf('GFP120X_%d_frame_%d.tif',mm,fstart-mm+1);
    GFP=imread([fdir1 fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',mm,((ii-1)*60)+mm); 
%             fnamer=sprintf('DAPI120X_%d_frame_%d.tif',mm,fstart-mm+1);
    RFP=imread([fdir1 fnamer]);
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    image=zeros(FLb,FLd,3);
    mycolors1=zeros(2001,3);
    mycolors1(1:2001,2)=0:0.0005:1;
    %     maxg=single(min(GFP(:))+1000);
    ming=single(min(GFP(:))+150);
    maxg=single(max(GFP(:)))-150; % -150 for 1 & 2
    %     maxg=2000; 
    minr=single(min(RFP(:)))+400;
    maxr=single(max(RFP(:)))-50; 
    GFPmask1=GFP>ming;
    GFPmask2=GFP<maxg;
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
    GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
    image(:,:,2)=GFP2;

    RFPmask1=RFP>minr;
    RFPmask2=RFP<maxr;
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
    RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
    mycolors2=zeros(2001,3);
    mycolors2(1:2001,1)=0:0.0005:1;
    image(:,:,1)=RFP2;
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    Image=rgb2gray(image);
    BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
    sz=size(BWn);
    c = kmeans(BWn(:), 2, 'MaxIter', 10000);
    CC=reshape(c,sz(1),sz(2));
    jk=mode(CC(:));
    BWsn=CC~=jk;
    %     BWsn = imbinarize(BWn,0.15);
    seD = strel('diamond',2);
    %BWfinal = imerode(BWnobord,seD);
    BWfinaln = imerode(BWsn,seD);
    BWfinaln = bwareaopen(BWfinaln, 200);
    Ln1 = imclearborder(BWfinaln, 4);
    windowSize = 51;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(Ln1), kernel, 'same');
    binaryImage = blurryImage > 0.5; 
    se90 = strel('line', 30, 90);
    se0 = strel('line', 30, 0);
    %     binaryImage = imdilate(binaryImage, [se90 se0]);
    %     Lnls=imresize(binaryImage,[356,458]);
    Lnls=imresize(BWfinaln,[QPIb-QPIa+1,QPId-QPIc+1]);
    BWNucl=zeros(512,512);
    BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls; % before FL magnifier
%             BWNucl(131:486,55:512)=Lnls; % after FL magnifier
    BWdilt=imdilate(BWNucl,[se0 se90]);
    BWCyto=BWfinal.*(1-BWNucl).*(BWdilt);
    figure(2);
    imoverlay(Abkg_mass(135:490,117:496,1),SGf(135:490,117:496).*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
    colormap(temp); hold on;
    [Bn,Lnlls] = bwboundaries(BWNucl(135:490,117:496),'noholes');
   for ee=1:length(Bn)
        boundary = Bn{ee};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
    end
%             DImage2=imresize(DImage,[356,458,3]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    figure(1)
    imshow(DImage3(135:490,117:496,:)); hold on;

    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=220; ybase=240; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textdown=sprintf('10 µm');
    text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
    pause(2);
end

