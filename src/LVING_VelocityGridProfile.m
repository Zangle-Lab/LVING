
%% LVING main code for computing the growth rate inside cells at pixel resolution
% SP Nov 2021
% first section: define image processing parameters (file locations, 
% processing options, etc.)
% second section: load phase images, the time captured
% third section: initialie grid parameters overlayed on image
% fourth section: background correct phase images to flaten
% fifth section: perform QPV using SSD on images, store displacements and
% deform initialized grid with the displacements
% sixth section: track the biomass change due to displacement of grid
% seventh section: compute growth rate at each CV using mass over time 
% computed in earlier section
% eigth section: store the results

clc; 
close all;
clear all;

%% set the parameters for computation
for mm=1:12   % set up for multiple fovs

    for rr=1:4     % for loop over multiple time points for each fov

        % define file path where phase data is stored
        fdir='S:\Data\Soorya\FixedCellImaging_May2018\MCF7_35mm_120X\Trial1\';

        %***SY****%
        tcg=1;                         % # of frames averaged over time
        xcg=4;                         % window size for spatial averaging of SSD, also size of CV for tracking
        gs=15;                         % SSD window size, chose an odd number
        sgap=1;                        % difference bwn frames undergoing SSD, must be less than (numf-tcg)
        % bgap=4;                      % difference between frames when computing velocity for grid interpolation correction
        numf=30;                       % total number of time points considered for calculation (interval GFP was imaged)
        medfiltsz=4*xcg;               % median filter size
        pxl_conv=0.21;                 % pixel size for 120X mag in um
        Ref_inc=0.18;                  % refractive index increment in um3/pg
        oplf=0.623;                    % wavelength of illumination light in um
        Pixel_area=pxl_conv*pxl_conv;  % pixel size in um
        w=3*gs;                        % SSD search window size
        %***SY***%

        fstart=((numf/2)*(rr-1))+1;

        %% stored the phase images and the time

        fname=sprintf('QPM120X_%d_frame_%d.mat',mm,fstart);      %first image which is position reference for GC calculation
        fname1=sprintf('QPM120X_%d_frame_%d.tif',mm,fstart);     % corresponding tif raw file used for time point metadata
        D=dir([fdir fname1]);
        Timenum(1)=D.datenum;
        lnumf=ceil(numf./sgap)-1;
        load([fdir fname],'Phase');
        for uu=1:ceil(numf./sgap)
            fname=sprintf('QPM120X_%d_frame_%d.mat',mm,(((uu-1)*sgap)+fstart));
            load([fdir fname],'Phase');
            fname1=sprintf('QPM120X_%d_frame_%d.tif',mm,(((uu-1)*sgap)+fstart));
            D=dir([fdir fname1]);
            Timenum(uu+1)=D.datenum;
            Time(uu)=(Timenum(uu+1)-Timenum(1))*24*60;   %in min
            D_stored2(:,:,uu)=Phase;
        end
        Tmass=reshape(sum(sum(D_stored2,1),2),1,ceil(numf./sgap));
        plot(Tmass)
        clear Abkg Phase L R B M;

        %% initialise the computation by defining the zero vectors to which the results are stored
        
        sz = size(D_stored2);
        kkf = sz(3)-(tcg);       % index of last frame to be used as CurrD image in velocity computation
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
        dX=zeros((sz(1)-gs),(sz(2)-gs),'single');
        dY=dX;
        clear X0 Y0 sX0 sY0; 
        
        %% compute mass of whole cell over time points considered

        A_stored2=zeros(sz,'single');
        Abkg_stored2=A_stored2;
        for ss=1:lnumf+1
            [B,A_stored2(:,:,ss)]=imagebackground_poly4(D_stored2(:,:,ss));
            Abkg_stored2(:,:,ss)= A_stored2(:,:,ss)*(oplf*Pixel_area/Ref_inc); 
            mass(ss)=sum(sum(Abkg_stored2(:,:,ss).*B));
        end
        
        %% compute the intracellular velocity and store results in terms of sptial grid interpolation over time
        
        for kk=1:kkf 
                CurrD = gpuArray(Abkg_stored2(:,:,kk));
                NextD = gpuArray(Abkg_stored2(:,:,kk+1));
                XCs=SSD_corr(NextD,CurrD,gs); 
                clear CurrD NextD;
                XCs=gather(XCs);
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
                    [xm,ym]=findvalley_v3(XCsp(:,:,ii,jj));
                    dX1(ii+2,jj+2) = xm;
                    dY1(ii+2,jj+2) = ym;
                end
            end
            clear XCsp; 
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

        %% use grid deformation to track biomass movement in space and time

        massThresh=0.001;   % mass threshold below which it is considered background
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
            Abkg_mass(:,:,kk+1)=CVtracking(xex,yey,xcg,Abkg_stored3);
            MassAm(kk+1)=sum(sum(Abkg_mass(:,:,kk+1)));
            fprintf('Completed CV loop %d\n',kk);
        end
        toc;

        %% fit the mass over time and find the slope of the growth curve to find the growth rate at each pixel

        for pp=1:sz(1)-(2*xcg) %(sz(1)-gs+1)/xcg
            for qq=1:sz(2)-(2*xcg) %(sz(2)-gs+1)/xcg
                xData=Time(1:kkf+1); %reshape(Abkg_massavg(pp,qq,:),1,kkf);
                yData=reshape(Abkg_mass(pp,qq,1:kkf+1),1,kkf+1);
                [fitData,coeff,const] = polyfit(xData, yData, 1);
                Massavg(pp,qq)=mean(yData);
                GC(pp,qq)=fitData(1)/const(2);
                if median(yData)~=0
                    GrowthData(pp,qq)=GC(pp,qq)/Massavg(pp,qq);
                else
                    GrowthData(pp,qq)=0;
                end
            end
        end
        imagesc(GC)

        %% validate and store results at mat files

        Mask1=Abkg_mass(:,:,1)>0.02;
        Mask2=GC>-1.2;
        [Rr,cft,cnst]=polyfit(Time',mass,1);
        GrCnt=(Rr(1)/cnst(2))/mean(mass);
        Gt=sum(sum(Mask1.*GC))./mean(MassAm);
        fname=sprintf('WS%d_cell%d.mat',rr,mm);
        save([fdir1 fname]);
        dXold=dX; dYold=dY;
        clearvars -except mm dXold dYold;

    end
end
