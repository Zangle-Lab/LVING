%% using farmework from Dr.Zangle code and computation from code 'Velocity_grid_profile_7'
%% rev 3 on 2 May 2018: SSD, pixel resolution velocity, spatial average of SSD, Find valley at minimum point, Control volume movement tracking
%% rev 5 on 31 May 2018: SSD calculation to reduce code time, vectorize SSD calculations, subtract whole images
%% rev 6 on 11 June 2018: no interpolation of velocity results, larger control volume of tracking
%% mistake in tracking rectified at line 99 and 100, grid interplotaion with tracking wrong in previous revisions
%% rev 13 introduced deforming grid on mass tracking instead of assumption of constant volume tracking
%% rev 13_5 on 7 June 2019: the alignment of interpolated grid with mass tracking has been corrected
%% rev 13_6 on 13th June 2019: Overlapping windows for control volume tracking
% clc; 
% close all;
% clear all;
% tic;
% function [GC]=VelocityGridProfile_rev13_4(rr)

%% set the parameters for computation
for mm=1:12
% mm=3,6;
%start again from 20, to 32
for rr=1:4
%     fdir1='Q:\Data\Soorya\U2OSFUCCIImaging120X_2019\P5_21Aug2019\MassGenResults_rev13_6\cell2\Results\';
%     fdir1='F:\Data\Soorya\MCF7Imaging_120X_May2019\P14_17May2019\MassGenTrackingResults_rev136_Trial1\cell7\demo\';
    fdir='S:\Data\Soorya\FixedCellImaging_May2018\MCF7_35mm_120X\Trial1\';

    %***SY****%
    tcg=1;    % # of frames averaged over time
    xcg=4;    % window size for spatial averaging of SSD, also size of CV for tracking
    gs=15;   % SSD window size, chose an odd number
    % bgs=35;
    sgap=1;   %difference bwn frames undergoing SSD, must be less than (numf-tcg)
    % bgap=4;   %difference between frames when computing velocity for grid interpolation correction
    numf=30;   % total number of frames considered for calculationGFP
    %     tgap=(30/60); %in minutes
    medfiltsz=4*xcg;
    pxl_conv=0.21;  % for 120X
    Ref_inc=0.18;    %in um3/pg
    oplf=0.623;
    Pixel_area=pxl_conv*pxl_conv;
    w=3*gs;    % SSD search window size
    %***SY***%

%%
    fstart=((numf/2)*(rr-1))+1;

    %% stored the phase images and the time
    %***SY***% 
%     Time=(tgap*(1:sgap:(numf)))./60;       % in hours, time corresponding to the phase image arrays
    fname=sprintf('QPM120X_%d_frame_%d.mat',mm,fstart);                 %first image which is position reference for GC calculation
    fname1=sprintf('QPM120X_%d_frame_%d.tif',mm,fstart);  
    D=dir([fdir fname1]);
    Timenum(1)=D.datenum;
    %***SY***%
    lnumf=ceil(numf./sgap)-1;
    load([fdir fname],'Phase');
%     [B,~]=imagebackground_poly4(Phase);
    for uu=1:ceil(numf./sgap)
        fname=sprintf('QPM120X_%d_frame_%d.mat',mm,(((uu-1)*sgap)+fstart));  %***SY***%
        load([fdir fname],'Phase');
        fname1=sprintf('QPM120X_%d_frame_%d.tif',mm,(((uu-1)*sgap)+fstart));
        D=dir([fdir fname1]);
        Timenum(uu+1)=D.datenum;
        Time(uu)=(Timenum(uu+1)-Timenum(1))*24*60;   %in min
%         D_stored2(:,:,uu)= imagebackground_poly4_kmeans(Abkg);
        D_stored2(:,:,uu)=Phase;
%         A_stored2(:,:,uu)=Abkg;
    %     D_stored2(:,:,uu)=-Phase-B;     % store background corrected phase images as an array
%         [L,R]=imagesegment_aggressive(Abkg);
%         Outline(:,:,uu)=L;            % draw outline of cell in which calculation is performed
    end
    Tmass=reshape(sum(sum(D_stored2,1),2),1,ceil(numf./sgap));
    plot(Tmass)
    clear Abkg Phase L R B M;

    % for qs=1:11 
    %     if Tmass(qs)<1800
    %         fname=sprintf('cell_4_frame_%d.mat',(((qs-1)*sgap)+fstart));  %***SY***%
    %         load([fdir fname],'Phase','Abkg');
    %         D_stored2(:,:,qs)= imagebackground_poly4_kmeans(Phase);
    %         Tmass(qs)=sum(sum(D_stored2(:,:,qs)));
    %     end 
    % end
    % D_stored2(:,:,10)=D_stored2(:,:,9);
    % Ds = load([fdir, 'data_allframes'], 'Time', 'D_stored2'); %load stored data workspace from plotHighResData.m
    % Ds.D_stored2 = single(Ds.D_stored2(:,:,1:100));
    % Ds.Time=Ds.Time(1:100,1);

    %% initialise the computation by defining the zero vectors to which the results are stored
    sz = size(D_stored2);
    kkf = sz(3)-(tcg);       % index of last frame to be used as CurrD image in velocity computation
%     bkf=floor(kkf/bgap);
    [X,Y] = meshgrid(1:sz(2), 1:sz(1));
    X = single(X);
    Y = single(Y);
    [X0,Y0]=meshgrid(1:sz(1)+1,1:sz(1)+1);
    sX0 = single(X0);
    sY0 = single(Y0);
    sX = single(sX0);
    sY = single(sY0);
    XS = single(sX0);
    YS = single(sY0);
%     bX = sX;
%     bY = sY;
%     dX1=zeros(floor((sz(1)-gs+1)/xcg),'single');
%     dY1=zeros(floor((sz(2)-gs+1)/xcg),'single');
    dX=zeros((sz(1)-gs),(sz(2)-gs),'single');
    dY=dX;
%     bdX=zeros(floor((sz(1)-bgs+1)/xcg),floor((sz(2)-bgs+1)/xcg),'single');
%     bdY=bdX;
    clear X0 Y0 sX0 sY0; 
    
%%
    A_stored2=zeros(sz,'single');
    Abkg_stored2=A_stored2;
%     Outline=A_stored2;
    for ss=1:lnumf+1
    %     [L,R]=imagesegment_aggressive(D_stored2(:,:,ss));
        [B,A_stored2(:,:,ss)]=imagebackground_poly4(D_stored2(:,:,ss));
        Abkg_stored2(:,:,ss)= A_stored2(:,:,ss)*(oplf*Pixel_area/Ref_inc);        
%         [Abkg_stored2(:,:,ss),BWmap]= imagebackground_poly4_kmeans(A_stored2(:,:,ss)*(oplf*Pixel_area/Ref_inc));
        mass(ss)=sum(sum(Abkg_stored2(:,:,ss).*B));
    end
%     Abkg_stored2=A_stored2.*(oplf*Pixel_area/Ref_inc);

%     for yy=(2*xcg)+1-(xcg/2):xcg:sz(1)-(2*xcg)+(xcg/2)
%         for xx=(2*xcg)+1-(xcg/2):xcg:sz(2)-(2*xcg)+(xcg/2)
%             Abkg_mass(((yy+1)/xcg)-1,((xx+1)/xcg)-1,1)=sum(sum(Abkg_stored2(yy:yy+xcg-1,xx:xx+xcg-1,1)));
%             Abkg_outline(((yy+1)/xcg)-1,((xx+1)/xcg)-1,1)=sum(sum(Outline(yy:yy+xcg-1,xx:xx+xcg-1,1)));
%         end
%     end

    %% Perform SSD on the small gap images and store the results
%     bk=1;
%     if rr>1
%         dX(:,:,1:(numf/2)-1)=dXold(:,:,(numf/2)+1:numf-1);
%         dY(:,:,1:(numf/2)-1)=dYold(:,:,(numf/2)+1:numf-1);
%         for kk=1:(numf/2)-1
%             for pp=medfiltsz+1:sz(1)-medfiltsz+1
%                 for qq=medfiltsz+1:sz(2)-medfiltsz+1
%                     sx=sX(pp,qq);
%                     rx=floor(sx);
%                     sy=sY(pp,qq);
%                     ry=floor(sy);
%                     if rx>0 && rx<sz(1)-medfiltsz
%                         if ry>0 && ry<sz(1)-medfiltsz
%                             sX(pp,qq) = sX(pp,qq) + (dX(ry,rx,kk));
%                             sY(pp,qq) = sY(pp,qq) + (dY(ry,rx,kk));
%                         end
%                     end
%                 end
%             end
%             XS(:,:,kk+1) = sX;
%             YS(:,:,kk+1) = sY;
%             fprintf('Completed loop %d\n',kk);
%         end
%         kstart=(kkf+1)/2;
%     else 
%         kstart=1;
%     end
    %%
    for kk=1:kkf 
%         if rem(kk,bgap)~=0
            CurrD = gpuArray(Abkg_stored2(:,:,kk));
            NextD = gpuArray(Abkg_stored2(:,:,kk+1));
%             CurrD=gpuArray(imtophat(single(D_stored2(:,:,kk)),strel('sphere',3)));
%             NextD=gpuArray(imtophat(single(D_stored2(:,:,kk+1)),strel('sphere',3)));
            XCs=SSD_corr_rev9(NextD,CurrD,gs); clear CurrD NextD;
            XCs=gather(XCs);
            XCsp=zeros((2*gs)+1,(2*gs)+1,sz(1),sz(2),'single');
            for uu=1:sz(1)-xcg+1    %spatial averaging of SSD
                for vv=1:sz(2)-xcg+1
                    XCsp(:,:,uu,vv)=mean(mean(XCs(:,:,uu:uu+xcg-1,vv:vv+xcg-1),3),4);
                end
            end
%           XCsp1=movmean(movmean(XCs,xcg,3),xcg,4);
%           XCsp=XCsp1(:,:,xcg-1:xcg:sz(1),xcg-1:xcg:sz(2));
            clear XCs;
%             dX=zeros(floor(sz(1)/xcg),floor(sz(2)/xcg),'single');
%             dY=dX;
%             wXCsp = libpointer('singlePtr',XCsp);
%             parfor ii =1:floor(sz(1)/xcg)
% %                 [dX1(ii,:),dY1(ii,:)]=FindValley(sz,xcg,wXCsp,ii);
%                 [dX1(ii,:),dY1(ii,:)]=FindValley(sz,xcg,squeeze(XCsp(:,:,ii,:)));
%             end
%           tic;
          dX1=zeros(sz(1)+1);
          dY1=dX1;
          for ii =1:sz(1)-3
               for jj =1:sz(2)-3
                  [xm,ym]=findvalley_v3_rev4(XCsp(:,:,ii,jj));
                  dX1(ii+2,jj+2) = xm;
                  dY1(ii+2,jj+2) = ym;
               end
          end
%         toc;
%         for uu=1:xcg:sz(1)-xcg+1    %spatial averaging of SSD
%             for vv=1:xcg:sz(2)-xcg+1
%                 dX2(((uu-1)/xcg)+1,((vv-1)/xcg)+1)=mean(mean(dX1(uu:uu+xcg-1,vv:vv+xcg-1),1),2);
%                 dY2(((uu-1)/xcg)+1,((vv-1)/xcg)+1)=mean(mean(dY1(uu:uu+xcg-1,vv:vv+xcg-1),1),2);
%             end
%         end
        clear XCsp; 
%         dX=zeros(sz(1)+1);
%         dY=dX;
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
        for pp=medfiltsz+1:sz(1)-medfiltsz+1
            for qq=medfiltsz+1:sz(2)-medfiltsz+1
                sx=sX(pp,qq);
                rx=floor(sx);
                sy=sY(pp,qq);
                ry=floor(sy);
                if rx>0 && rx<sz(1)-medfiltsz
                    if ry>0 && ry<sz(1)-medfiltsz
                        sX(pp,qq) = sX(pp,qq) + (dX(ry,rx,kk));
                        sY(pp,qq) = sY(pp,qq) + (dY(ry,rx,kk));
                    end
                end
            end
        end
        XS(:,:,kk+1) = sX;
        YS(:,:,kk+1) = sY;
        fprintf('Completed loop %d\n',kk);
        clear CurrD NextD dX1 dY1;
    end
%     else
%         CurrD = gpuArray(imtophat(single(D_stored2(:,:,bk)),strel('sphere',3)));
%         NextD = gpuArray(imtophat(single(D_stored2(:,:,bk+bgap)),strel('sphere',3)));
%         XCs=SSD_corr_rev9(NextD,CurrD,bgs);
%         XCs=gather(XCs);
%         XCsp=zeros((2*bgs)+1,(2*bgs)+1,sz(1)/xcg,sz(2)/xcg,'single');
%         for uu=1:xcg:sz(1)-xcg+1    %spatial averaging of SSD
%             for vv=1:xcg:sz(2)-xcg+1
%                 XCsp(:,:,((uu-1)/xcg)+1,((vv-1)/xcg)+1)=mean(mean(XCs(:,:,uu:uu+xcg-1,vv:vv+xcg-1),3),4);
%             end
%         end
%         clear XCs;
%         bdX1=zeros(floor(sz(1)/xcg),floor(sz(2)/xcg),'single');
%         bdY1=bdX1;
% %         wXCsp = libpointer('singlePtr',XCsp);
% %         parfor ii =1:floor(sz(1)/xcg)
% % %             [bdX1(ii,:),bdY1(ii,:)]=FindValley(sz,xcg,wXCsp.Value);
% %             [bdX1(ii,:),bdY1(ii,:)]=FindValley(sz,xcg,squeeze(XCsp(:,:,ii,:)));
% %         end
%         for ii =1:sz(1)/xcg
%             for jj =1:sz(2)/xcg
%                 [xm,ym]=findvalley_v3_rev4(XCsp(:,:,ii,jj));
%                 bdX1(ii,jj) = xm;
%                 bdY1(ii,jj) = ym;
%             end
%         end
%         clear XCsp;
%         for ee=medfiltsz+1:floor(sz(1)/xcg)-medfiltsz
%             for ff=medfiltsz+1:floor(sz(2)/xcg)-medfiltsz
%                 bUx=bdX1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
%                 bux=median(bUx(:),'omitnan');
%                 bUy=bdY1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
%                 buy=median(bUy(:),'omitnan');
%                 if abs(bdX1(ee,ff))>(1.2*bux) || abs(bdX1(ee,ff))<(0.8*bux)
%                     bdX(ee,ff,((bk-1)/bgap)+1)=bux;
%                 else
%                     bdX(ee,ff,((bk-1)/bgap)+1)=bdX1(ee,ff);
%                 end
%                 if abs(bdY1(ee,ff))>(1.2*buy) || abs(bdY1(ee,ff))<(0.8*buy)
%                     bdY(ee,ff,((bk-1)/bgap)+1)=buy;
%                 else
%                     bdY(ee,ff,((bk-1)/bgap)+1)=bdY1(ee,ff);
%                 end
%             end
%         end
%         for pp=1:floor(sz(1)/xcg)
%             for qq=1:floor(sz(1)/xcg)
%                 bx=bX(pp,qq);
%                 rx=floor(bx);
%                 by=bY(pp,qq);
%                 ry=floor(by);
%                 if rx>0 && rx<(sz(1)/xcg)-xcg
%                     if ry>0 && ry<(sz(1)/xcg)-xcg
%                        bX(pp,qq) = bX(pp,qq) + (bdX(ry,rx,((bk-1)/bgap)+1)./xcg);
%                        bY(pp,qq) = bY(pp,qq) + (bdY(ry,rx,((bk-1)/bgap)+1)./xcg);
%                     end
%                 end
%             end
%         end
%         XS(:,:,kk+1) = sX;
%         YS(:,:,kk+1) = sY;
%         sX=bX;
%         sY=bY;
%         fprintf('Completed loop %d\n',kk);
%         clear bdX1 bdY1;
%         bk=kk+1;
%         end

    % display the results, comment if not required
%     MB = imresize(Outline(9:504,9:504,1),[124 124]);
%     DX=(1-MB).*dX(1:124,1:124,1);
%     DY=(1-MB).*dY(1:124,1:124,1);
%     % % for uu=1:4:sz(1)-4+1
%     % %     for vv=1:4:sz(2)-4+1
%     % %         S(((uu-1)/4)+1,((vv-1)/4)+1)=mean2(DX(uu:uu+3,vv:vv+3));
%     % %         Q(((uu-1)/4)+1,((vv-1)/4)+1)=mean2(DY(uu:uu+3,vv:vv+3));
%     % %     end
%     % % end
%     x=[(-sz(2)/2) (sz(2)/2)];                                           %x axis ticks with origin as center
%     y=[(-sz(1)/2) (sz(1)/2)]; 
%     [Vxo,Vyo]=meshgrid((-sz(2)/xcg)+xcg+1:(xcg/2):(sz(2)/xcg)-xcg,(-sz(1)/xcg)+xcg+1:(xcg/2):(sz(1)/xcg)-xcg);   % center axis at origin to quiver
%     imagesc(x,y,Abkg_stored2(:,:,1));
%     hold on 
%     %quiver(Vxo,Vyo,dX(1:124,1:124,1),dY(1:124,1:124,1),3);
%     quiver(Vxo,Vyo,dX(:,:,1),dY(:,:,1),3);

%% 
tic;
    massThresh=0.001;   % mass below which it is considered background, ***SY***%
%     Abkg_OL=Abkg_outline>0;
%     tic; for pp=(2*xcg)+1-(xcg/2):sz(1)-xcg-(xcg/2)
%         for qq=(2*xcg)+1-(xcg/2):sz(2)-xcg-(xcg/2)
%             if Abkg_stored2(pp,qq,1)<massThresh
%                 Abkg_stored3(pp-xcg-(xcg/2),qq-xcg-(xcg/2),1)=0;
%             else
%                 Abkg_stored3(pp-xcg-(xcg/2),qq-xcg-(xcg/2),1)=Abkg_stored2(pp,qq,1);
%             end
%         end
%     end; toc;
    Abkg_mass=zeros(sz(1)-(2*xcg),sz(2)-(2*xcg),kkf+1,'single');
    MassAm=zeros(kkf+1,1,'single');
    mass=MassAm;
    Abkg_mask2=Abkg_stored2(:,:,1)>massThresh;
    Abkg_stored3=(Abkg_stored2(:,:,1)).*Abkg_mask2;
    SOPL1=movsum(movsum(Abkg_stored3,xcg,1),xcg,2);
    Abkg_mass(:,:,1)=(SOPL1(xcg+1:sz(1)-xcg,xcg+1:sz(2)-xcg))/(xcg*xcg);
    % Abkg_massch=zeros(sz(1));
    % Abkg_massavg=zeros(sz(1));
    mass(1)=sum(sum(Abkg_stored3(xcg+1:sz(1)-xcg,xcg+1:sz(2)-xcg,1)));
    MassAm(1)=sum(sum(Abkg_mass(:,:,1)));
    for kk=1:kkf
        Abkg_mask2=Abkg_stored2(:,:,kk+1)>massThresh;
        Abkg_stored3=(Abkg_stored2(:,:,kk+1)).*Abkg_mask2;
        mass(kk+1)=sum(sum(Abkg_stored3));  
        SOPL=Abkg_stored3;
        xex=XS(:,:,kk+1);
        yey=YS(:,:,kk+1);
        Abkg_mass(:,:,kk+1)=CVtracking_rev1(xex,yey,xcg,Abkg_stored3);
        MassAm(kk+1)=sum(sum(Abkg_mass(:,:,kk+1)));
        fprintf('Completed CV loop %d\n',kk);
    end
    toc;

    %% fitting the mass track and finding the slope of the growth curve to find the growth rate at each pixel

   for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
        for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
            xData=Time(1:kkf+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
            yData=reshape(Abkg_mass(pp,qq,1:kkf+1),1,kkf+1);
            [fitData,coeff,const] = polyfit(xData, yData, 1);
            Massavg(pp,qq)=mean(yData);
%             GC(pp,qq)=fitData1(1);
            GC(pp,qq)=fitData(1)/const(2);
            if median(yData)~=0
                GrowthData(pp,qq)=GC(pp,qq)/Massavg(pp,qq);
            else
                GrowthData(pp,qq)=0;
            end
        end
    end
    % GCfull=interp2(sX0,sY0,GC,X,Y);
    imagesc(GC)
    % time_elapsed=toc;

%     %% actual growth rate to compare
%     figure(2);
%     plot(Time,mass);
%     for zz=1:kkf+1
%         Outline1(:,:,zz)=imagesegment_kmeans(Abkg_stored2(:,:,zz));
%         Outadj(:,:,zz)=imresize(Outline1(7:506,7:506,zz),[125,125]);
%         Abkg_massadj(:,:,zz)=Abkg_mass(:,:,zz).*Outadj(:,:,zz);
%         Abkg_storedadj(:,:,zz)=Abkg_stored2(:,:,zz).*Outline1(:,:,zz);
%         Massadj(zz)=sum(sum(Abkg_storedadj(:,:,zz)));
%     end
%     [Rr,cft,cnst]=polyfit(Time,mass,1);
%     MassAvg=mean(mass);
%     GrCnt=(Rr(1)/cnst(2))/MassAvg;   % for the whole cell
%     Growth=(sum(sum(GrowthData.*Abkg_mass(:,:,1))))./(sum(sum(Abkg_mass(:,:,1))));% from tracking
%     Gt=sum(sum(Outadj(:,:,1).*GC))/sum(sum(mean(Abkg_massadj,3),2),1);
%     fname=sprintf('WS%d_cell1.mat',rr);
%     save([fdir1 fname]);
%     clear all     n 
%% validation
%     for tt=1:kkf+1 
        Mask1=Abkg_mass(:,:,1)>0.02;
%         Mask1(200:250,270:300)=1;
        Mask2=GC>-1.2;
%         MassAm(tt)=sum(sum(Mask1.*Mask2.*Abkg_mass(:,:,tt)));
%     end
    [Rr,cft,cnst]=polyfit(Time',mass,1);
    GrCnt=(Rr(1)/cnst(2))/mean(mass);
%     Gt=sum(sum(Mask1.*Mask2.*GC))/mean(MassAm);
    Gt=sum(sum(Mask1.*GC))./mean(MassAm);
    fname=sprintf('WS%d_cell%d.mat',rr,mm);
    save([fdir1 fname]);
    dXold=dX; dYold=dY;
    clearvars -except mm dXold dYold;
end
end

%% compare growth rate constant from mass tracking inside a cell and the overall mass change in the cell
% for oo=1:240 
%     B=imagebackground_poly4(D_stored2(:,:,oo));
%     Abkg_mass2(:,:,oo)=D_stored2(:,:,oo)-B;
% end
% Growth=(sum(sum(GC(10:502,10:502).*Abkg_mass(10:502,10:502,1))))./(sum(sum(Abkg_mass(10:502,10:502,1))));
% Time=(0.5:0.5:120)./60;       % in hours, time corresponding to the phase image arrays
% mass=sum(sum(Abkg_mass2(:,:,1)));
% for ii=2:240
%     mass(ii)=sum(sum(Abkg_mass2(:,:,ii)));
% %     mass_ch(ii)=(mass(ii)-mass(ii-1))./(Time(ii)-Time(ii-1));
% %     massavg(ii)=(mass(ii)+mass(ii-1))./2;
% end
% coeff=polyfit(Time,mass,1);
% totalMassAvg=mean(mass);
% GrowthWhole=coeff/totalMassAvg;