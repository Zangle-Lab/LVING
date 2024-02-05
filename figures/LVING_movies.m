
%% LVING movies 
clc; clear all;

%% video M1: RPE cell cycle LVING growth map over time

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
%fdirFL='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
fdirFL = 'S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\Trial2\';
%fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdirM='S:\Data\Soorya\RPEFUCCIImaging_2020\Pn2_31Jan2020\MassGenResults_rev136\Results_2hr\';
cellno =2;
for dd=7:30
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM]);
    DD = Abkg_stored2(:,:,1);
%     [BWfinal,L] = BWmask_initial(PreMaskImage);
%     [BWfinalF,L1] = BWmask_initial(Abkg_stored2(:,:,120));
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',cellno,60*(dd-1)+cellno);
    GFP=imread([fdir fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',cellno,60*(dd-1)+cellno);
    RFP=imread([fdir fnamer]);
    
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    
    limg=[single(min(min(GFP)))+250,single(max(max(GFP)))];
    limr=[single(min(min(RFP)))+150,single(max(max(RFP)))];
    GFPmask1=GFP>limg(1);
    GFPmask2=GFP<limg(2);
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limg(2))+((1-GFPmask1).*limg(1));
    GFP2=(GFP1-limg(1))/(limg(2)-limg(1));
    RFPmask1=RFP>limr(1);
    RFPmask2=RFP<limr(2);
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*limr(2))+((1-RFPmask1).*limr(1));
    RFP2=(RFP1-limr(1))/(limr(2)-limr(1));

    SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    
    figure(1);
    subplot(1,2,1);
    imoverlay(Abkg_mass(:,:,1),SGf,[-0.0005, 0.0005],[],parula, 0.2, gca); 
    colormap(map); %ylim([59 339]); xlim([160 440]);
    
    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=443; ybase=490; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textup=sprintf('%d hr',(dd-6));
    text('units','pixels','position',[20 30],'fontsize',20,'color','w','string',textup)
    textdown=sprintf('10 µm');
    text('units','pixels','position',[410 40],'fontsize',20,'color','w','string',textdown)
    hold off
    
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    subplot(1,2,2);
    imshow(DImage3); hold on; % ylim([210 490]); xlim([150 430]);
    pause(1);
    
    set(gcf,'color','w'); pause(1);
    fdir='T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Videos\M1_RPECellCycle\';
    Fname=sprintf('RPECellCycle_frame%d.jpg',dd);
    M=getframe(gcf);
    imwrite(M.cdata, [fdir Fname])
    
end

%% video M2: MCF7 cell cycle LVING growth map over time
len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
%fdirFL='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
fdirFL = 'F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\Trial1\';
%fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdirM='F:\Data\Soorya\MCF7Imaging_120X_May2019\P15_20May2019\MassGenResults_rev136\Results_2hr\';
cellno =3;
for dd=13:22
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM]);
    DD = Abkg_stored2(:,:,1);
%     [BWfinal,L] = BWmask_initial(PreMaskImage);
%     [BWfinalF,L1] = BWmask_initial(Abkg_stored2(:,:,120));
    fnameg=sprintf('GFP120X_%d_frame_%d.tif',cellno,fstart);
    GFP=imread([fdir fnameg]); 
    fnamer=sprintf('RFP120X_%d_frame_%d.tif',cellno,fstart);
    RFP=imread([fdir fnamer]);
    
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    
    limg=[single(min(min(GFP)))+40,single(max(max(GFP)))-50];
    limr=[single(min(min(RFP)))+10,single(max(max(RFP)))-60];
    GFPmask1=GFP>limg(1);
    GFPmask2=GFP<limg(2);
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limg(2))+((1-GFPmask1).*limg(1));
    GFP2=(GFP1-limg(1))/(limg(2)-limg(1));
    RFPmask1=RFP>limr(1);
    RFPmask2=RFP<limr(2);
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*limr(2))+((1-RFPmask1).*limr(1));
    RFP2=(RFP1-limr(1))/(limr(2)-limr(1));

    SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    
    figure(1);
    subplot(1,2,1);
    imoverlay(Abkg_mass(116:385,140:419,1),SGf(116:385,140:419),[-0.003, 0.003],[],parula, 0.2, gca); 
    colormap(map); %ylim([59 339]); xlim([160 440]);
    
    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=233; ybase=260; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textup=sprintf('%d hr',(dd-8));
    text('units','pixels','position',[15 20],'fontsize',20,'color','w','string',textup)
    textdown=sprintf('10 µm');
    text('units','pixels','position',[380 40],'fontsize',20,'color','w','string',textdown)
    hold off
    
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    subplot(1,2,2);
    imshow(DImage3(112:389,136:423,:)); hold on; % ylim([210 490]); xlim([150 430]);
    pause(1);
    
    set(gcf,'color','w'); pause(1);
    fdir='T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Videos\M2_MCF7CellCycle\';
    Fname=sprintf('MCFCellCycle_frame%d.jpg',dd);
    M=getframe(gcf);
    imwrite(M.cdata, [fdir Fname])
    
end

%% video M3: RPE cell autophagy, growth inhibition over time

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

QPIa=131; QPIb=486; QPIc=55; QPId=512; % after FL magnifier
FLa=1; FLb=512; FLc=1; FLd=671;
%fdirFL='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
fdirFL = 'Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\Trial1\';
%fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdirM='Q:\Data\Soorya\RPEAutophagyImaging_2020\Pn8_5March2020\MassGenResults_rev136\Results_2hr\';
cellno = 2;
for dd=2:30
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM]);
    DD = Abkg_stored2(:,:,1);
%     [BWfinal,L] = BWmask_initial(PreMaskImage);
%     [BWfinalF,L1] = BWmask_initial(Abkg_stored2(:,:,120));
    fnameg=sprintf('mAG120X_%d_frame_%d.tif',cellno,fstart);
    GFP=imread([fdirFL fnameg]); 
    fnamer=sprintf('mKO2120X_%d_frame_%d.tif',cellno,fstart);
    RFP=imread([fdirFL fnamer]);
    
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    
    limg=[single(min(min(GFP)))+250,single(max(max(GFP)))];
    limr=[single(min(min(RFP)))+150,single(max(max(RFP)))];
    GFPmask1=GFP>limg(1);
    GFPmask2=GFP<limg(2);
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limg(2))+((1-GFPmask1).*limg(1));
    GFP2=(GFP1-limg(1))/(limg(2)-limg(1));
    RFPmask1=RFP>limr(1);
    RFPmask2=RFP<limr(2);
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*limr(2))+((1-RFPmask1).*limr(1));
    RFP2=(RFP1-limr(1))/(limr(2)-limr(1));

    SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    
    figure(1);
    subplot(1,2,1);
    imoverlay(Abkg_mass(:,:,1),SGf,[-0.0005, 0.0005],[],parula, 0.2, gca); 
    colormap(map); %ylim([59 339]); xlim([160 440]);
    
    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=443; ybase=490; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textup=sprintf('%d hr',(dd-2));
    text('units','pixels','position',[20 30],'fontsize',20,'color','w','string',textup)
    textdown=sprintf('10 µm');
    text('units','pixels','position',[410 40],'fontsize',20,'color','w','string',textdown)
    hold off
    
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,1)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    subplot(1,2,2);
    imshow(DImage3); hold on; % ylim([210 490]); xlim([150 430]);
    pause(1);
    
    set(gcf,'color','w'); pause(1);
    fdir='T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Videos\M3_RPEAutophagy\';
    Fname=sprintf('RPEAutophagy_frame%d.jpg',dd);
    M=getframe(gcf);
    imwrite(M.cdata, [fdir Fname])
    
end

%% video M4: MCF7 cell autophagy growth inhibition over time

len = size(get(gcf, 'Colormap'), 1);
map=temp(len);

QPIa=112; QPIb=389; QPIc=136; QPId=495;  % before FL magnifier
FLa=1; FLb=512; FLc=1; FLd=688;
%fdirFL='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\Trial1\';
fdirFL = 'Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\Trial1\';
%fdirM='R:\Data\Soorya\RPEFUCC1_2020\Pn7_2Sep2020\MassGenResults_rev136\Results_2hr\';
fdirM='Q:\Data\Soorya\MCF7AutophagyImaging120X_2019\P17_20Dec2019\MassGenResults_rev136\Results_2hr\';
cellno = 5;
for dd=2:14
    fnameM=sprintf('WS%d_cell%d.mat',dd, cellno); load([fdirM fnameM]);
    DD = Abkg_stored2(:,:,1);
%     [BWfinal,L] = BWmask_initial(PreMaskImage);
%     [BWfinalF,L1] = BWmask_initial(Abkg_stored2(:,:,120));
    fnameg=sprintf('GFP120X_%d_frame_%d.tif',cellno,fstart);
    GFP=imread([fdir fnameg]); 
    fnamer=sprintf('DAPI120X_%d_frame_%d.tif',cellno,fstart);
    RFP=imread([fdir fnamer]);
    
    RFP=RFP(:,1:FLd); GFP=GFP(:,1:FLd);
    
    limg=[single(min(min(GFP)))+10,single(max(max(GFP)))-100];
    limr=[single(min(min(RFP)))+10,single(max(max(RFP)))-100];
    GFPmask1=GFP>limg(1);
    GFPmask2=GFP<limg(2);
    GFP1=(single(GFP).*GFPmask1.*GFPmask2)+((1-GFPmask2).*limg(2))+((1-GFPmask1).*limg(1));
    GFP2=(GFP1-limg(1))/(limg(2)-limg(1));
    RFPmask1=RFP>limr(1);
    RFPmask2=RFP<limr(2);
    RFP1=(single(RFP).*RFPmask1.*RFPmask2)+((1-RFPmask2).*limr(2))+((1-RFPmask1).*limr(1));
    RFP2=(RFP1-limr(1))/(limr(2)-limr(1));

    SGf=zeros(512,512); SGf = imfilter((GC), fspecial('gaussian', [50 50], 1));
    SM = Abkg_mass(:,:,1)>0.001;
    SGf = SGf.*SM;

    figure(1);
    subplot(1,2,1);
    imoverlay(Abkg_mass(116:385,140:419,1),SGf(116:385,140:419),[-0.0005, 0.0005],[],parula, 0.2, gca); 
    colormap(map); %ylim([59 339]); xlim([160 440]);
    
    barsize=10; pxlsize=0.238/1000;
    bar_2 = barsize./pxlsize./1000./2; %half-width of scalebar, in pixels
    hold on; xbase=233; ybase=260; 
    H = plot([xbase-bar_2 xbase+bar_2], ybase+[0 0], '-w', 'LineWidth', 2);
    textup=sprintf('%d hr',(dd-2));
    text('units','pixels','position',[15 20],'fontsize',20,'color','w','string',textup)
    textdown=sprintf('10 µm');
    text('units','pixels','position',[380 40],'fontsize',20,'color','w','string',textdown)
    hold off
    
    DImage=zeros(QPIb-QPIa+1,QPId-QPIc+1,3); 
    DImage(:,:,3)=imresize(RFP2,[QPIb-QPIa+1,QPId-QPIc+1]); 
    DImage(:,:,2)=imresize(GFP2,[QPIb-QPIa+1,QPId-QPIc+1]);
    DImage3=zeros(512,512,3); DImage3(QPIa:QPIb,QPIc:QPId,:)=DImage;
    subplot(1,2,2);
    imshow(DImage3(112:389,136:423,:)); hold on; % ylim([210 490]); xlim([150 430]);
    pause(1);
    
    set(gcf,'color','w'); pause(1);
    fdir='T:\Data\Soorya\LIVINGPaperFigures_2021\Rev5_Videos\M4_MCF7Autophagy\';
    Fname=sprintf('MCF7Autophagy_frame%d.jpg',dd);
    M=getframe(gcf);
    imwrite(M.cdata, [fdir Fname])
    
end
