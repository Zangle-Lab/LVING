%% LVING setup figure elements
clear all; clc;

%% figure 1a : Quantitative phase image of an RPE cell, colormap temp

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\QPM120X_2_frame_782.mat');
[M,Abkg1,SS]=imagebackground_poly4_double(Phase*(0.623*(0.21*0.21)/0.18));
figure(1);
imagesc(Abkg1(127:497,77:462)./(0.21*0.21)); colormap(temp); axis off image; colorbar; 
caxis([-0.1 1]); set(gcf,'color','w');
c.FontSize=18; set(gcf,'color','w');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
textup=sprintf('0 min');
text('units','pixels','position',[20 320],'fontsize',15,'color','w','string',textup)

%% figure 1b : difference between QPI images of RPE cell over 10 minutes

load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\QPM120X_2_frame_792.mat');
[M,Abkg2,SS]=imagebackground_poly4_double(Phase*(0.623*(0.21*0.21)/0.18));
figure(2);
imagesc((Abkg2(127:497,77:462)-Abkg1(127:497,77:462))./(0.21*0.21));
colormap(temp); axis off image; colorbar; caxis([-0.4 0.4]);
set(gcf,'color','w');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-k', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','k','string',textdown)
textup=sprintf('10-0 min');
text('units','pixels','position',[20 320],'fontsize',15,'color','k','string',textup)

%% figure 1c : intracellular velocity map on the RPE cell

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\WS14_cell2.mat')

figure(3);
dX2=imresize(dX(128:497,78:462,3),[376/8,390/8]);
% dX2=imresize(dX(167:421,195:425),[254/8,230/8]); 
dY2=imresize(dY(128:497,78:462,3),[376/8,390/8]);
[Vxo,Vyo]=meshgrid(-178:8:206,-128:8:240);
x=[-178 206];
y=[-128 240];
[M,~]=imagebackground_poly4(D_stored2(:,:,3));
MB=imresize(M(128:497,78:462),0.125);
DX=MB.*dX2;
DY=MB.*dY2;
DX(31,2)=0.238*8;
Vavg=sqrt((dX(:,:,3).^2)+(dY(:,:,3).^2));
imagesc(x,y,Vavg(128:497,78:462)/(0.21*0.21)); colormap hsv; %(200:450,150:350)
axis off image; set(gcf,'color','w'); caxis([0 100]);
hold on;
quiver(Vxo,Vyo,DX,DY,3,'k');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=140; ybase=150; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[90 180],'fontsize',15,'color','w','string',textdown)
textup=sprintf('1 µm/min');
text('units','pixels','position',[200 90],'fontsize',15,'color','k','string',textup)

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWdilt=imdilate(M,[se0 se90]);
hold on;
[Bn,Lnlls] = bwboundaries(BWdilt(128:497,78:462),'noholes');
Bn{1,1}(:,1)=Bn{1,1}(:,1)-128; Bn{1,1}(:,2)=Bn{1,1}(:,2)-178;
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end

figure(4);
imagesc(x,y,Abkg_stored2(128:497,78:462,1)/(0.21*0.21)); colormap(temp); %(200:450,150:350)
axis off image; set(gcf,'color','w'); caxis([-0.1 1]);
hold on;
quiver(Vxo,Vyo,DX,DY,3,'k');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=140; ybase=150; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[90 180],'fontsize',15,'color','w','string',textdown)
textup=sprintf('1 µm/min');
text('units','pixels','position',[200 90],'fontsize',15,'color','k','string',textup)


%% figure 1d : selected CV tracks over 30 minutes overlayed on the RPE cell QPI image 

clear all;
load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\WS14_cell2.mat')
len = size(get(gcf, 'Colormap'), 1);
temp=temp(len);
imagesc(Abkg_stored2(:,:,1)/(0.21*0.21)); axis off image; colormap(temp); caxis([-0.1 1]);

zlength=30;
TT(1,1)=384; TT(1,2)=312;
TT(2,1)=359; TT(2,2)=280;
TT(3,1)=316; TT(3,2)=263;
TT(4,1)=312; TT(4,2)=319;

for rr=1:length(TT)
    tt=TT(rr,1); ww=TT(rr,2);
    hold on   
    xos(1,:)=reshape(XS(tt,ww,1:zlength),1,zlength); yos(1,:)=reshape(YS(tt,ww,1:zlength),1,zlength);
    xos(2,:)=reshape(XS(tt+4,ww,1:zlength),1,zlength); yos(2,:)=reshape(YS(tt+4,ww,1:zlength),1,zlength);
    xos(3,:)=reshape(XS(tt+4,ww+4,1:zlength),1,zlength); yos(3,:)=reshape(YS(tt+4,ww+4,1:zlength),1,zlength);
    xos(4,:)=reshape(XS(tt,ww+4,1:zlength),1,zlength); yos(4,:)=reshape(YS(tt,ww+4,1:zlength),1,zlength);

    for xx=1:4
        plot(xos(xx,:),yos(xx,:),'k', 'MarkerSize', 5, 'LineWidth', 1)
        pause(0.2); 
    end
    hold on;
    h1=rectangle('Position',[xos(1,1) yos(1,1) (xos(4,1)-xos(1,1)) (yos(3,1)-yos(1,1))],'EdgeColor','r'); h1.LineWidth=2;
    for Pn=10:10:30
        line([xos(1,Pn), xos(2,Pn)], [yos(1,Pn), yos(2,Pn)], 'Color', 'g', 'LineWidth', 1);
        line([xos(2,Pn), xos(3,Pn)], [yos(2,Pn), yos(3,Pn)], 'Color', 'g', 'LineWidth', 1);
        line([xos(3,Pn), xos(4,Pn)], [yos(3,Pn), yos(4,Pn)], 'Color', 'g', 'LineWidth', 1);
        line([xos(4,Pn), xos(1,Pn)], [yos(4,Pn), yos(1,Pn)], 'Color', 'g', 'LineWidth', 1);
    end
end
ylim([128 497]);
xlim([78 462]);

%% figure 1e : slope of CV mass vs time plot 

MassCV=zeros(4,zlength);
for RR=1:length(TT)
    pp=TT(RR,1); qq=TT(RR,2);
    MassCV(RR,:)=Abkg_mass(qq,pp,1:zlength);
    [p,S] = polyfit(Time(1:zlength),MassCV(RR,:),1);
    f = polyval(p,Time(1:zlength));
    plot(Time(1:zlength),MassCV(RR,:),'o',Time(1:zlength),f,'-');
    hold on;
end

%% figure 1f : computed LVING growth map (colormap parula) overlayed on the RPE cell QPI image

load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\WS14_cell2.mat')
SGf=zeros(512,512);
SGfin = imfilter((GC), fspecial('gaussian', [50 50], 1));
SGf(5:508,5:508)=SGfin;
SM = Abkg_stored2(:,:,1)>0.001;
SGf = SGf.*SM;
% SGf(SGf==0) = NaN;
figure(4);
imoverlay(Abkg_stored2(128:497,78:462,1),SGf(128:497,78:462).*60,[-0.06, 0.06],[],parula, 0.2, gca); 
colorbar; axis image; set(gcf,'color','w');
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=200; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)

%% figure 2a and 2b : G1 phase LVING growth map & fluorescence images

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
  
load('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\MassGenResults_rev136\Results_2hr\WS6_cell2.mat')
SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.01; SGf = SGf.*SM; SGf(SGf==0) = NaN;
figure(1);
imoverlay(Abkg_mass(:,:,1),SGf,[-0.001, 0.001],[],parula, 0.2, gca); 
colormap(map); ylim([59 339]); xlim([160 440]);
GFP=imread('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\Trial2\mAG120X_2_frame_302.tif');
RFP=imread('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\Trial2\mKO2120X_2_frame_302.tif');
minr=200; ming=200; maxr=4000; maxg=2000;
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
ylim([59 339]); xlim([160 440]);

%% figure 2c and 2d : G1-S phase transition LVING growth map & fluorescence images

load('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\MassGenResults_rev136\Results_2hr\WS26_cell1.mat')
SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.01; SGf = SGf.*SM; SGf(SGf==0) = NaN;
figure(3);
imoverlay(Abkg_mass(:,:,1),SGf,[-0.001, 0.001],[],parula, 0.2, gca); 
colormap(map); ylim([80 360]); xlim([120 400]);
GFP=imread('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\Trial2\mAG120X_1_frame_1501.tif');
RFP=imread('S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\Trial2\mKO2120X_1_frame_1501.tif');
minr=200; ming=200; maxr=4000; maxg=2000;
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
figure(4)
imshow(DImage3); hold on; ylim([80 360]); xlim([120 400]);

%% figure 2e and 2f : S phase LVING growth map & fluorescence images

load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\WS14_cell2.mat')
SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.01; SGf = SGf.*SM; SGf(SGf==0) = NaN;
figure(5);
imoverlay(Abkg_mass(:,:,1),SGf,[-0.001, 0.001],[],parula, 0.2, gca); 
colormap(map); ylim([210 490]); xlim([150 430]);
GFP=imread('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\mAG120X_2_frame_782.tif');
RFP=imread('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\mKO2120X_2_frame_782.tif');
minr=200; ming=200; maxr=4000; maxg=2000;
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
figure(6)
imshow(DImage3); hold on; ylim([210 490]); xlim([150 430]);

%% figure 2g and 2h : G2 phase LVING growth map & fluorescence images

load('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\WS18_cell2.mat')
SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
SM = Abkg_mass(:,:,1)>-0.01; SGf = SGf.*SM; SGf(SGf==0) = NaN;
figure(7);
imoverlay(Abkg_mass(:,:,1),SGf,[-0.001, 0.001],[],parula, 0.2, gca); 
colormap(map); ylim([180 460]); xlim([190 470]);
GFP=imread('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\mAG120X_2_frame_1022.tif');
RFP=imread('R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\mKO2120X_2_frame_1022.tif');
minr=200; ming=200; maxr=4000; maxg=2000;
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
figure(8)
imshow(DImage3); hold on; ylim([180 460]); xlim([190 470]);
  
%% figure 2i : specific growth in RPE cells in different cell cycle phases

sgrN0_all=zeros(180,4); sgrC0_all=zeros(180,4); sgrW0_all=zeros(180,4); CellPhase_all=zeros(180,4);

load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat')
sgrN0_all(1:30,:)=sgrN0.*InDepth; sgrC0_all(1:30,:)=sgrC0.*InDepth; sgrW0_all(1:30,:)=sgrW0.*InDepth; CellPhase_all(1:30,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Water2.3uL_22May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat')
sgrN0_all(31:60,:)=sgrN0.*InDepth; sgrC0_all(31:60,:)=sgrC0.*InDepth; sgrW0_all(31:60,:)=sgrW0.*InDepth; CellPhase_all(31:60,:)=CellPhase.*InDepth;

load('S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat')
sgrN0_all(61:90,:)=sgrN0.*InDepth; sgrC0_all(61:90,:)=sgrC0.*InDepth; sgrW0_all(61:90,:)=sgrW0.*InDepth; CellPhase_all(61:90,:)=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn10_Water2.3uL_18June2020\MassGenresults_rev136\Results_2hr\NCWsgr_rev6.mat')
sgrN0_all(91:120,:)=sgrN0.*InDepth; sgrC0_all(91:120,:)=sgrC0.*InDepth; sgrW0_all(91:120,:)=sgrW0.*InDepth; CellPhase_all(91:120,:)=CellPhase.*InDepth;
sgrN0_all(isinf(abs(sgrN0_all))|isnan(sgrN0_all))=0; sgrC0_all(isinf(abs(sgrC0_all))|isnan(sgrC0_all))=0; sgrW0_all(isinf(abs(sgrW0_all))|isnan(sgrW0_all))=0;

G1WS=(CellPhase_all==1); G1SWS=(CellPhase_all==2); SWS=(CellPhase_all==3); G2WS=(CellPhase_all==4);
SRNcG1_all=G1WS.*sgrN0_all; SRNcG1S_all=G1SWS.*sgrN0_all; SRNcS_all=SWS.*sgrN0_all; SRNcG2_all=G2WS.*sgrN0_all;
SRCyG1_all=G1WS.*sgrC0_all; SRCyG1S_all=G1SWS.*sgrC0_all; SRCyS_all=SWS.*sgrC0_all; SRCyG2_all=G2WS.*sgrC0_all;
SRWG1_all=G1WS.*sgrW0_all; SRWG1S_all=G1SWS.*sgrW0_all; SRWS_all=SWS.*sgrW0_all; SRWG2_all=G2WS.*sgrW0_all;

for cc=1:length(CellPhase_all)
    SRNcG1(cc,1)=mean(nonzeros(SRNcG1_all(cc,:))); 
    SRNcG1S(cc,1)=mean(nonzeros(SRNcG1S_all(cc,:)));
    SRNcS(cc,1)=mean(nonzeros(SRNcS_all(cc,:)));  
    SRNcG2(cc,1)=mean(nonzeros(SRNcG2_all(cc,:)));
    SRCyG1(cc,1)=mean(nonzeros(SRCyG1_all(cc,:))); 
    SRCyG1S(cc,1)=mean(nonzeros(SRCyG1S_all(cc,:)));
    SRCyS(cc,1)=mean(nonzeros(SRCyS_all(cc,:))); 
    SRCyG2(cc,1)=mean(nonzeros(SRCyG2_all(cc,:)));
    SRWG1(cc,1)=mean(nonzeros(SRWG1_all(cc,:)));
    SRWG1S(cc,1)=mean(nonzeros(SRWG1S_all(cc,:)));
    SRWS(cc,1)=mean(nonzeros(SRWS_all(cc,:))); 
    SRWG2(cc,1)=mean(nonzeros(SRWG2_all(cc,:)));
end
SRNcG1(isnan(SRNcG1))=0; SRNcG1S(isnan(SRNcG1S))=0; SRNcS(isnan(SRNcS))=0; SRNcG2(isnan(SRNcG2))=0;
SRCyG1(isnan(SRCyG1))=0; SRCyG1S(isnan(SRCyG1S))=0; SRCyS(isnan(SRCyS))=0; SRCyG2(isnan(SRCyG2))=0;
SRWG1(isnan(SRWG1))=0; SRWG1S(isnan(SRWG1S))=0; SRWS(isnan(SRWS))=0; SRWG2(isnan(SRWG2))=0;

NdifG1mn=median(nonzeros(SRNcG1)); NdifG1Smn=median(nonzeros(SRNcG1S)); NdifSmn=median(nonzeros(SRNcS)); NdifG2mn=median(nonzeros(SRNcG2));
NdifG1std=std(SRNcG1(SRNcG1~=0))/sqrt(nnz(SRNcG1)); NdifG1Sstd=std(SRNcG1S(SRNcG1S~=0))/sqrt(nnz(SRNcG1S)); NdifSstd=std(SRNcS(SRNcS~=0))/sqrt(nnz(SRNcS)); NdifG2std=std(SRNcG2(SRNcG2~=0))/sqrt(nnz(SRNcG2));
CdifG1mn=median(nonzeros(SRCyG1)); CdifG1Smn=median(nonzeros(SRCyG1S)); CdifSmn=median(nonzeros(SRCyS)); CdifG2mn=median(nonzeros(SRCyG2));
CdifG1std=std(SRCyG1(SRCyG1~=0))/sqrt(nnz(SRCyG1)); CdifG1Sstd=std(SRCyG1S(SRCyG1S~=0))/sqrt(nnz(SRCyG1S)); CdifSstd=std(SRCyS(SRCyS~=0))/sqrt(nnz(SRCyS)); CdifG2std=std(SRCyG2(SRCyG2~=0))/sqrt(nnz(SRCyG2));
WdifG1mn=median(nonzeros(SRWG1)); WdifG1Smn=median(nonzeros(SRWG1S)); WdifSmn=median(nonzeros(SRWS)); WdifG2mn=median(nonzeros(SRWG2));
WdifG1std=std(SRWG1(SRWG1~=0))/sqrt(nnz(SRWG1)); WdifG1Sstd=std(SRWS(SRWG1S~=0))/sqrt(nnz(SRWG1S)); WdifSstd=std(SRWS(SRWS~=0))/sqrt(nnz(SRWS)); WdifG2std=std(SRWG2(SRWG2~=0))/sqrt(nnz(SRWG2));

Ny=[WdifG1mn,WdifG1Smn,WdifSmn,WdifG2mn;NdifG1mn,NdifG1Smn,NdifSmn,NdifG2mn;CdifG1mn,CdifG1Smn,CdifSmn,CdifG2mn];
% Ny=[NdifG1mn,NdifSmn,NdifG2mn;CdifG1mn,CdifSmn,CdifG2mn;WdifG1mn,WdifSmn,WdifG2mn];
Nerr=[WdifG1std,WdifG1Sstd,WdifSstd,WdifG2std;NdifG1std,NdifG1Sstd,NdifSstd,NdifG2std;CdifG1std,CdifG1Sstd,CdifSstd,CdifG2std];
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
ydata{1,4}=nonzeros(SRWG2); ydata{2,4}=nonzeros(SRNcG2); ydata{3,4}=nonzeros(SRCyG2); 

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

set(gcf,'Color','w'); ylim([-0.3 0.3]);
ylabel('Specific growth rate (/hr)');

%% figure 2j and 2k: specific growth in nucl & cyto and GFP-RFP FL intensity
% along cell cycle growth time

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
fdirFL='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';

cellno =2;
for dd=2:16
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM]);
    PreMaskImage = Abkg_stored2(:,:,1);
    [BWfinal,L] = BWmask_initial(PreMaskImage);
    [BWfinalF,L1] = BWmask_initial(Abkg_stored2(:,:,120));
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',cellno,60*(dd-1)+cellno);
    GFP=imread([fdir fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',cellno,60*(dd-1)+cellno);
    RFP=imread([fdir fnamer]);
    limg=[single(min(min(GFP)))+150,single(max(max(GFP)))];
    limr=[single(min(min(RFP)))+150,single(max(max(RFP)))];
    GFPmask1=GFP>limg(1);
    GFPmask2=GFP<limg(2);
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limg(2))+((1-GFPmask1).*limg(1));
    GFP2=(GFP1-limg(1))/(limg(2)-limg(1));
    RFPmask1=RFP>limr(1);
    RFPmask2=RFP<limr(2);
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*limr(2))+((1-RFPmask1).*limr(1));
    RFP2=(RFP1-limr(1))/(limr(2)-limr(1));
    Image=zeros(512,688,3); Image(:,:,1)=RFP2; Image(:,:,2)=GFP2;
    Image=rgb2gray(Image);
    BWn = imfilter(Image, fspecial('gaussian', [10 10], 1));
    szF=size(BWn);
    c = kmeans(BWn(:), 2, 'MaxIter', 10000);
    CC=reshape(c,szF(1),szF(2));
    jk=mode(CC(:));
    BWsn=CC~=jk;
    seD = strel('diamond',2);
    BWfinaln = imerode(BWsn,seD);
    BWfinaln = bwareaopen(BWfinaln, 200);
    Ln1 = imclearborder(BWfinaln, 4);
    windowSize = 51;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(Ln1), kernel, 'same');
    binaryImage = blurryImage > 0.5; 
    se90 = strel('line', 8, 90);
    se0 = strel('line', 8, 0);
    binaryImage = imdilate(binaryImage, [se90 se0]);
    %     Lnls=imresize(binaryImage,[356,458]);
    Lnls=imresize(binaryImage,[QPIb-QPIa+1,QPId-QPIc+1]);
    Lnls=imclearborder(Lnls,4);
    se90 = strel('line', 12, 90);
    se0 = strel('line', 12, 0);
    BWfinalC=imerode(BWfinal,[se90 se0]);
    BWNucl=zeros(512,512);
    BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls;
    BWNucl=BWNucl.*BWfinalC;
    BWCyto=((BWfinalC-BWNucl)==1);
    SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    BWNC=abs(SGf)>0.0003;
%     imoverlay(Abkg_mass(:,:,1),SGf,[-0.001, 0.001],[],parula, 0.2, gca); 
%     pause();
    MassW=sum(sum(BWfinalC.*Abkg_stored2(:,:,1)));
    MassN=sum(sum(BWNucl.*Abkg_stored2(:,:,1)));
    MassC=sum(sum(BWCyto.*Abkg_stored2(:,:,1)));
    %BWsum=sum(sum(BWfinalF));
    sgrN0(dd)=((sum(sum(SGf.*BWNucl(5:508,5:508).*BWNC)))./MassN)*60;
    sgrC0(dd)=((sum(sum(SGf.*BWCyto(5:508,5:508).*BWNC)))./MassC)*60;
    sgrW0(dd)=((sum(sum(SGf.*BWfinalC(5:508,5:508).*BWNC)))./MassW)*60;
end
plot(sgrN0); hold on; 
plot(sgrW0); plot(sgrC0);


fdir='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
WSno =35;
cellno=2;
% fdir ='K:\Data\Soorya\RPEFUCCIImaging_2020\Pn7_28Feb2020\Trial1\';
% WSno=100;
% cellno=3;
for dd=1:WSno
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',cellno,30*(dd-1)+cellno);
    GFP=imread([fdir fnameg]); 
    %figure(1); imagesc(GFP); pause();
    IntG(dd)=sum(sum(GFP));
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',cellno,30*(dd-1)+cellno);
    RFP=imread([fdir fnamer]);
    %figure(2); imagesc(RFP); pause();
    IntR(dd)=sum(sum(RFP));
end
plot(IntG/IntG(1)); hold on; plot(IntR/IntR(WSno));

fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
for SS=12:45
    fnameM=sprintf('WS%d_cell7.mat',SS);
    load([fdirM fnameM]);
    imagesc(Abkg_stored2(:,:,1)); pause(1);
end

%% figure 3a: puncta overlay image of RPE cell, G1 phase

clear all; clc;

fdir='S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\MassGenResults_rev136\Results_2hr\';
fdir1='S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\Trial1\';
mm=12; ii=2;

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
imoverlay(Abkg_mass(135:430,117:416,1),SGf(135:430,117:416).*60,[-0.08, 0.08],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWdilt(135:430,117:416),'noholes');
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

%% figure 3b: puncta overlay image of RPE cell, S phase.
clear all; clc;
fdir='S:\Data\Soorya\RPEFUCCIImaging_2020\Pn6_29Aug2020\MassGenResults_rev136\Results_2hr\';
fdir1='S:\Data\Soorya\RPEFUCCIImaging_2020\Pn6_29Aug2020\Trial1\';
mm=2; ii=30;

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
imoverlay(Abkg_mass(:,:,1),SGf.*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
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


%% figure 4: a and b: Drug treatment dox results
clear all; clc;

fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn3Doxorubicin_21March2020\Trial1\';
mm=3; ii=3;

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
figure(1);
imoverlay(Abkg_mass(135:430,117:416,1),SGf(135:430,117:416).*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWNucl(135:430,117:416),'noholes');
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
figure(2)
imshow(DImage3(135:430,117:416,:)); hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)

% figure 4 d and e: Drug treatment homo results
clear all; clc;

fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.36um_3May2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.36um_3May2020\Trial1\';
mm=13; ii=3;

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
figure(1);
imoverlay(Abkg_mass(135:430,117:416,1),SGf(135:430,117:416).*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
colormap(temp); hold on;
[Bn,Lnlls] = bwboundaries(BWNucl(135:430,117:416),'noholes');
for ee=1:length(Bn)
    boundary = Bn{ee};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
figure(2)
imshow(DImage3(135:430,117:416,:)); hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)


%% figure 3 panel 4: Cytoplasm growth inhibiting drug Homoharringtonine treatment overall results

load('K:\Data\Soorya\RPEDrugTest_2020\Pn6_DMSO2.1uL_4May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1=sgrN0.*InDepth(:,1:3); sgrN0_1(1,4)=0; sgrC0_1=sgrC0.*InDepth(:,1:3); sgrC0_1(1,4)=0; sgrW0_1=sgrW0.*InDepth(:,1:3); sgrW0_1(1,4)=0; CellPhase_1=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn10_DMSO2.1uL_18June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1(31:60,:)=sgrN0.*InDepth; sgrC0_1(31:60,:)=sgrC0.*InDepth; sgrW0_1(31:60,:)=sgrW0.*InDepth; CellPhase_1(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.036um_3May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2=sgrN0.*InDepth; sgrC0_2=sgrC0.*InDepth; sgrW0_2=sgrW0.*InDepth; CellPhase_2=CellPhase.*InDepth;
load('T:\Data\Soorya\RPEDrugTest2020\Pn7_Homoharringtonine0.036uM_13Jun2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2(31:60,:)=sgrN0.*InDepth; sgrC0_2(31:60,:)=sgrC0.*InDepth; sgrW0_2(31:60,:)=sgrW0.*InDepth; CellPhase_2(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn6_Homoharringtonine0.36um_3May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3=sgrN0.*InDepth; sgrC0_3=sgrC0.*InDepth; sgrW0_3=sgrW0.*InDepth; CellPhase_3=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn8_Homoharringtonine0.36uM_14Jun2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3(31:60,:)=sgrN0.*InDepth; sgrC0_3(31:60,:)=sgrC0.*InDepth; sgrW0_3(31:60,:)=sgrW0.*InDepth; CellPhase_3(31:60,:)=CellPhase.*InDepth;

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
set(gca, 'XTickLabel',{'DMSO','Homo 0.036','Homo 0.36'});
ylim([-0.15 0.15]);


%% figure S13 f,g,h: Nuclear sgr in G1, S, G2 phases 

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
set(gca, 'XTickLabel', {'DMSO','Homo 0.036','Homo 0.36'});


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
set(gca, 'XTickLabel',{'DMSO','Homo 0.036','Homo 0.36'});

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
set(gca, 'XTickLabel', {'DMSO','Homo 0.036','Homo 0.36'});

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
set(gca, 'XTickLabel', {'DMSO','Homo 0.036','Homo 0.36'});

%% figure S15 : Cycloheximide drug impact on RPE cells

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

% figure a,b,c,d
% for cyclo 0.01
fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.1um_21May2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.1um_21May2020\Trial1\';
mm=20; ii=4;

% for cyclo 0.53
fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide0.53um_20May2020\MassGenResults_rev136\Results_2hr\';
fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide0.53um_20May2020\Trial1\';
mm=7; ii=4;


% for cyclo 2.88
% fdir='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide2.88um_20May2020\MassGenResults_rev136\Results_2hr\';
% fdir1='K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide2.88um_20May2020\Trial1\';
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
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    figure(1)
    imshow(DImage3(135:430,117:416,:)); hold on;

    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=220; ybase=240; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textdown=sprintf('10 µm');
    text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
    pause(2);
end

% figure S15 e, f, g, h

load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1=sgrN0.*InDepth; sgrC0_1=sgrC0.*InDepth; sgrW0_1=sgrW0.*InDepth;  CellPhase_1=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn10_Ethanol0.04uL_18June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_1(31:60,:)=sgrN0.*InDepth; sgrC0_1(31:60,:)=sgrC0.*InDepth; sgrW0_1(31:60,:)=sgrW0.*InDepth; CellPhase_1(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.1um_21May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2=sgrN0.*InDepth; sgrC0_2=sgrC0.*InDepth; sgrW0_2=sgrW0.*InDepth; CellPhase_2=CellPhase.*InDepth;
load('T:\Data\Soorya\RPEDrugTest2020\Pn7_Cycloheximide0.1uM_11June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_2(31:60,:)=sgrN0.*InDepth; sgrC0_2(31:60,:)=sgrC0.*InDepth; sgrW0_2(31:60,:)=sgrW0.*InDepth; CellPhase_2(31:60,:)=CellPhase.*InDepth;
load('K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide0.53um_20May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3=sgrN0.*InDepth; sgrC0_3=sgrC0.*InDepth; sgrW0_3=sgrW0.*InDepth; CellPhase_3=CellPhase.*InDepth;
load('S:\Data\Soorya\RPEDrugTest_2020\Pn7_Cycloheximide0.53uM_11June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
sgrN0_3(31:60,:)=sgrN0.*InDepth; sgrC0_3(31:60,:)=sgrC0.*InDepth; sgrW0_3(31:60,:)=sgrW0.*InDepth; CellPhase_3(31:60,:)=CellPhase.*InDepth;
% load('K:\Data\Soorya\RPEDrugTest_2020\Pn6_Cycloheximide2.88um_20May2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
% sgrN0_4=sgrN0.*InDepth; sgrN0_4(1,4)=0; sgrC0_4=sgrC0.*InDepth; sgrC0_4(1,4)=0; sgrW0_4=sgrW0.*InDepth; sgrW0_4(1,4)=0; CellPhase_4=CellPhase.*InDepth;
% load('T:\Data\Soorya\RPEDrugTest2020\Pn7_Cycloheximide2.88uM_11June2020\MassGenResults_rev136\Results_2hr\NCWsgr_rev6.mat');
% sgrN0_4(31:60,:)=sgrN0.*InDepth; sgrC0_4(31:60,:)=sgrC0; sgrW0_4(31:60,:)=sgrW0; CellPhase_4(31:60,:)=CellPhase;

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
set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53'});
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
set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53'});


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
set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53'});

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
set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53'});

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
set(gca, 'XTickLabel', {'Ethanol','Cyclo 0.01','Cyclo 0.53'});

%% figure 5 :  Autophagy of RPE cells

% figure 5 panel a to c: Autophagy growth map over time in RPE
clear all; clc;
len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
fdirRes='Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\MassGenResults_rev136\Results_2hr\';
fdirFL='Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\';
FrameEE=[2,10,20]; mm=2; %2, 6,10,11,27
QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;

for zz=1:length(FrameEE)
    fname=sprintf('WS%d_cell%d.mat',FrameEE(zz),mm);
    load([fdirRes fname],'GC','sz','fstart','Abkg_mass');
    SGf=zeros(512,512);
    SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    % SGf = S.GC;
    SM = Abkg_mass(:,:,1)>-0.01;
    SGf = SGf.*SM;
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',mm,fstart);
%             fnameg=sprintf('GFP120X_%d_frame_%d.tif',mm,fstart-mm+1);
    GFP=imread([fdirFL fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',mm,fstart); 
%             fnamer=sprintf('DAPI120X_%d_frame_%d.tif',mm,fstart-mm+1);
    RFP=imread([fdirFL fnamer]);
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    image=zeros(FLb,FLd,3);
    mycolors1=zeros(2001,3);
    mycolors1(1:2001,2)=0:0.0005:1;
    %     maxg=single(min(GFP(:))+1000);
    ming=single(min(GFP(:))+50);
    maxg=single(max(GFP(:)))-100; 
    %     maxg=2000; 
    minr=single(min(RFP(:)))+50;
    maxr=single(max(RFP(:)))-100; 
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
    BWNucl=zeros(512,512); DImageZ=zeros(512,512,3);
    BWNucl(QPIa:QPIb,QPIc:QPId)=Lnls;
    DImageZ(QPIa:QPIb,QPIc:QPId,:)=DImage;
%     figure(1); imshow(DImageZ(135:430,117:416,:));
% %     imoverlay(Abkg_mass(:,:,1),SGf,[-0.0005, 0.0005],[],parula, 0.2, gca); %colormap(map);
%     barsize=10; pxlsize=0.238/1000;
%     bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
%     hold on; xbase=220; ybase=240; 
%     H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
%     textdown=sprintf('10 µm');
%     text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
%     axis image;
    figure(zz);
    imoverlay(Abkg_mass(135:430,117:416,1),SGf(135:430,117:416).*60,[-0.03, 0.03],[],parula, 0.2, gca); %colormap(map);
    colormap(temp); hold on;
%     [Bn,Lnlls] = bwboundaries(BWNucl(135:430,117:416),'noholes');
%     hold on
%     for ee=1:length(Bn)
%         boundary = Bn{ee};
%         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
%     end
end

hold on;
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)


% figure % d : growth structure deacy rate
QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
fdirRes='Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\MassGenResults_rev136\Results_2hr\';
fdirFL='Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\';
cellno=2;
FrameEE=1:20;
for zz=1:length(FrameEE)
    fname=sprintf('WS%d_cell%d.mat',FrameEE(zz),cellno);
    load([fdirRes fname],'GC','sz','fstart','Abkg_mass','Abkg_stored2');
    DD=Abkg_stored2(:,:,1);
    SSf = imfilter(abs(DD), fspecial('gaussian', [20 20], 2));
    [junk threshold] = edge(SSf, 'sobel');
    fudgeFactor = 0.5; %was 0.4 for RBCs
    BWs = edge(SSf,'sobel', threshold * fudgeFactor);
    se90 = strel('line', 8, 90);
    se0 = strel('line', 8, 0);
    BWsdil = imdilate(BWs, [se90 se0]);
    BWdfill = imfill(BWsdil, 'holes');
    seD = strel('diamond',5);
    BWfinal = imerode(BWdfill,seD);
    BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
    se90 = strel('line', 5, 90);
    se0 = strel('line', 5, 0);
    BWfinal=imdilate(BWfinal,[se0 se90]);
    BWfinal=imclearborder(BWfinal,4);
    SGf=zeros(512,512);
    SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    normImage = mat2gray(SGf);
    level = graythresh(normImage);  % use in case the threshold goes wrong due to normalization with high negative pixels in cell background
    level=0.27; % level of threshold kept same as puncta degrades and there will be nothing left to threshold
    %BW = imbinarize(normImage,level-0.05);
    BW = normImage>=level;
    PunctaMag(zz) = (sum(sum(BW.*SGf.*BWfinal(5:508,5:508)))/sum(sum(DD.*BWfinal)))*60;
end

% figure 5. e to f: Overlap of LC3 puncta with degradtion puncta
% earlier figure(rev1) was cell 10
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
load('Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\QPM120X_5_frame_61.mat');
GFP=imread('Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\GFP120X_5_frame_61.tif');
DAPI=imread('Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\DAPI120X_5_frame_61.tif');
DAPI=DAPI(:,1:FLd); GFP=GFP(:,1:FLd);
image=zeros(FLb,FLd,3);
mycolors1=zeros(2001,3);
mycolors1(1:2001,2)=0:0.0005:1;
%     maxg=single(min(GFP(:))+1000);
ming =220; %ming=single(min(GFP(:))+50);
maxg =260; %maxg=single(max(GFP(:)))-150; % -150 for 1 & 2
%     maxg=2000; 
minr = 130; %minr=single(min(DAPI(:)))+30;
maxr =200; %maxr=single(max(DAPI(:)))-200; 
GFPmask1=GFP>ming;
GFPmask2=GFP<maxg;
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
image(:,:,2)=GFP2;

RFPmask1=DAPI>minr;
RFPmask2=DAPI<maxr;
RFP1=(single(DAPI).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,3)=RFP2;
DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
DImage(:,:,3)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
imshow(DImage);
[~,Abkg]=imagebackground_poly4(Phase);
DD=imfilter(Abkg, fspecial('gaussian', [20 20], 1));
SSf = imfilter(abs(DD), fspecial('gaussian', [20 20], 2));
[junk threshold] = edge(SSf, 'sobel');
fudgeFactor = 0.3; %was 0.4 for RBCs
BWs = edge(SSf,'sobel', threshold * fudgeFactor);
se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
seD = strel('diamond',5);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);
figure(1)
imshow(DImage); hold on;
[B,Llls] = bwboundaries(BWfinal(QPIa:QPIb,QPIc:QPId),'noholes');
hold on
for ee=1:length(B)
    boundary = B{ee};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
len = size(get(gcf, 'Colormap'), 1);
map=temp(len);
figure(2); imagesc(Abkg(QPIa:QPIb,QPIc:QPId).*(0.623/0.18)); colormap(temp); axis off image; colorbar; caxis([-0.1 1.5]); 
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
figure(3); 
load('Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\MassGenResults_rev136\Results_2hr\WS13_cell5.mat');

SGf=zeros(512,512);
SGfin = imfilter((GC), fspecial('gaussian', [50 50], 1));
SGf(5:508,5:508)=SGfin;
SM = Abkg_stored2(:,:,1)>0.001;
SGf = SGf.*SM;
% SGf(SGf==0) = NaN;
imoverlay(Abkg_stored2(QPIa:QPIb,QPIc:QPId,1),SGf(QPIa:QPIb,QPIc:QPId).*60,[-0.02, 0.02],[],parula, 0.2, gca);% colormap(map);
colorbar; axis image; set(gcf,'color','w'); colormap(temp);

% figure 5.h : Timing for degradation puncta and LC3 puncta
fdirM ='Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\MassGenResults_rev136\Results_2hr\';
fdirFL = 'Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\';
cellno=5;
frameEE=15;
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
for zz=1:frameEE
    fnameWS = sprintf('WS%d_cell%d.mat',zz,cellno);
    load([fdirM fnameWS]);
    fnameGFP = sprintf('GFP120X_%d_frame_%d.tif',cellno,fstart);
    GFP=imread([fdirFL fnameGFP]);
    GFP=GFP(:,1:FLd);
    
    DD=imfilter(Abkg_stored2(:,:,1), fspecial('gaussian', [20 20], 1));
    SSf = imfilter(abs(DD), fspecial('gaussian', [20 20], 2));
    [junk threshold] = edge(SSf, 'sobel');
    fudgeFactor = 0.3; %was 0.4 for RBCs
    BWs = edge(SSf,'sobel', threshold * fudgeFactor);
    se90 = strel('line', 8, 90);
    se0 = strel('line', 8, 0);
    BWsdil = imdilate(BWs, [se90 se0]);
    BWdfill = imfill(BWsdil, 'holes');
    seD = strel('diamond',5);
    BWfinal = imerode(BWdfill,seD);
    BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
    se90 = strel('line', 5, 90);
    se0 = strel('line', 5, 0);
    BWfinal=imdilate(BWfinal,[se0 se90]);
    
    %GFP2 = imfilter((GFP), fspecial('gaussian', [50 50], 1));
    ming =220; %ming=single(min(GFP(:))+50);
    maxg =250; %maxg=single(max(GFP(:)))-150; % -150 for 1 & 2 
    GFPmask1=GFP>ming;
    GFPmask2=GFP<maxg;
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
    GFP2=(GFP1-ming)/(maxg-ming);
    GFP2a =imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    GFP3 = imfilter((GFP2a), fspecial('gaussian', [50 50], 1));
    normImage = mat2gray(GFP3);
    level = graythresh(normImage);  % use in case the threshold goes wrong due to normalization with high negative pixels in cell background
    %level=0.55;
    BWLC3 = imbinarize(normImage,level+0.05);
    LC3Area(zz) = sum(sum(BWLC3.*BWfinal(QPIa:QPIb,QPIc:QPId)));
    
    levelGC= -0.001;
    BWGC = GC<=levelGC;
    GCMag(zz)=sum(sum(BWGC.*GC));
end
plot(LC3Area./LC3Area(1)); hold on; plot(GCMag./GCMag(1));

%% figure S16 : Western blot and densitometry analysis

%% figure S17: Mass vs time LC2 staining

% figure S17 a for RPE cell
Mass_all=zeros(60,1439);
load('Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_4March2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat','MassWhole','TimeWhole');
Mass_all(1:30,1:1155)=MassWhole;
load('Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat','MassWhole','TimeWhole');
Mass_all(31:60,1:1439)=MassWhole;
Mass_all(Mass_all==0) = NaN;
figure(2);
[lineOut, fillOut] = stdshade(Mass_all,0.5,'b',TimeWhole(1,:)./60,5);
ylim([150 400]); set(gcf,'color','w'); 
xlabel('Time (hr)'); ylabel('Mass (pg)');

% figure S17 RPE LC2 stain image
QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
load('Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\QPM120X_1_frame_631.mat');
GFP=imread('Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\mAG120X_1_frame_631.tif');
DAPI=imread('Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\mKO2120X_1_frame_631.tif');
DAPI=DAPI(:,1:FLd); GFP=GFP(:,1:FLd);
image=zeros(FLb,FLd,3);
mycolors1=zeros(2001,3);
mycolors1(1:2001,2)=0:0.0005:1;
%     maxg=single(min(GFP(:))+1000);
ming=single(min(GFP(:))+60);
maxg=single(max(GFP(:)))-70; % -150 for 1 & 2
%     maxg=2000; 
minr=single(min(DAPI(:)))+30;
maxr=single(max(DAPI(:)))-100; 
GFPmask1=GFP>ming;
GFPmask2=GFP<maxg;
GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*maxg)+((1-GFPmask1).*ming);
GFP2=(GFP1-ming)/(maxg-ming); % normalize GFP image
image(:,:,2)=GFP2;

RFPmask1=DAPI>minr;
RFPmask2=DAPI<maxr;
RFP1=(single(DAPI).*RFPmask1.*RFPmask2)+((1-RFPmask2).*maxr)+((1-RFPmask1).*minr);
RFP2=(RFP1-minr)/(maxr-minr); % normalize RFP image
mycolors2=zeros(2001,3);
mycolors2(1:2001,1)=0:0.0005:1;
image(:,:,1)=RFP2;
DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
imshow(DImage);
[~,Abkg]=imagebackground_poly4(Phase);
DD=imfilter(Abkg, fspecial('gaussian', [2 2], 1));
SSf = imfilter(abs(DD), fspecial('gaussian', [20 20], 2));
[junk threshold] = edge(SSf, 'sobel');
fudgeFactor = 0.6; %was 0.4 for RBCs
BWs = edge(SSf,'sobel', threshold * fudgeFactor);
se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
seD = strel('diamond',5);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);
figure(1)
imshow(DImage); hold on;
[B,Llls] = bwboundaries(BWfinal(QPIa:QPIb,QPIc:QPId),'noholes');
hold on
for ee=1:length(B)
    boundary = B{ee};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); %draw cytoplasmic boundary on growth data
end
barsize=10; pxlsize=0.238/1000;
bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
hold on; xbase=220; ybase=240; 
H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
textdown=sprintf('10 µm');
text('units','pixels','position',[240 35],'fontsize',15,'color','w','string',textdown)
pause(2);

% figure S17 C for MCF7 cells
Mass_all=zeros(46,1439);
load('Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat','MassWhole','TimeWhole');
Mass_all(1:15,1:1340)=MassWhole(1:15,:);
load('T:\Data\Soorya\MCF7AutophagyImaging120X_2019\P15_21Jan2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat','MassWhole','TimeWhole');
Mass_all(16:23,1:1439)=MassWhole(1:8,:);
load('T:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_21Dec2019\MassGenResult_rev136\Results_2hr\CellMassvsTime_rev1.mat','MassWhole','TimeWhole');
Mass_all(24:38,1:1253)=MassWhole(1:15,:);
load('T:\Data\Soorya\MCF7AutophagyImaging120X_2019\P15_21Jan2020\MassGenResults_rev136\Results_2hr\CellMassvsTime_rev1.mat')
Mass_all(39:46,1:1439)=MassWhole(1:8,:);
Mass_all(Mass_all==0) = NaN;
[lineOut, fillOut] = stdshade(Mass_all,0.5,'b',TimeWhole(1,:)./60,5);
ylim([150 400]); xlim([0 20]); set(gcf,'color','w'); 
xlabel('Time (hr)'); ylabel('Mass (pg)');

% figure S17 d : MCF7 Lc2 stain intensity over time
fdir='T:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_21Dec2019\Trial1\';
PunctaDensity=zeros(15,42); 
QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
for cellno=1:15
    CellPunctaDen=zeros(1,42);
    for frame=1:30:1231
       fnameGFP=sprintf('GFP120X_%d_frame_%d.tif',cellno,frame);
       GFP=imread([fdir fnameGFP]);
       fnameQ=sprintf('QPM120X_%d_frame_%d.mat',cellno,frame);
       load([fdir fnameQ]);
       [BWfinal,SS] = imagebackground_poly4(Phase);
       limG=single([min(min(GFP))+50,max(max(GFP))]);
%        GFP0 = imfilter(GFP, fspecial('gaussian', [1 1], 1));
%        GFPmask1=GFP0>limG(1);
%        GFPmask2=GFP0<limG(2);
%        GFP1=(single(GFP0).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limG(2))+((1-GFPmask1).*limG(1));
       GFP2=(single(GFP)-limG(1))./(limG(2)-limG(1)); % normalize GFP image
%        GFP3=imtophat(GFP2, strel('sphere',15));
       GFP3=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
       CellPunctaDen(1,((frame-1)/30)+1)=mean2(nonzeros(GFP3.*BWfinal(QPIa:QPIb,QPIc:QPId)));
    end
    PunctaDensity(cellno,:)=CellPunctaDen./CellPunctaDen(1,1);
end

fdir='Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\';
for cellno=1:15
    CellPunctaDen=zeros(1,42);
    for frame=1:30:1231
       fnameGFP=sprintf('GFP120X_%d_frame_%d.tif',cellno,frame);
       GFP=imread([fdir fnameGFP]);
       limG=single([min(min(GFP))+50,max(max(GFP))]);
       GFP0 = imfilter(GFP, fspecial('gaussian', [1 1], 1));
       GFPmask1=GFP0>limG(1);
       GFPmask2=GFP0<limG(2);
       GFP1=(single(GFP0).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limG(2))+((1-GFPmask1).*limG(1));
       GFP2=(single(GFP1)-limG(1))./(limG(2)-limG(1)); % normalize GFP image
%        GFP3=imtophat(GFP2, strel('sphere',15));
       CellPunctaDen(1,((frame-1)/30)+1)=mean2(nonzeros(GFP2));
    end
    PunctaDensity(cellno+15,:)=CellPunctaDen./CellPunctaDen(1,1);
end

[lineOut, fillOut] = stdshade(PunctaDensity,0.5,'g',0:0.5:20.5,5);
xlim([0 20]); ylim([0.5 2.5]); set(gcf,'color','w'); 
xlabel('Time (hr)'); ylabel('Lysossome FL intensity');

