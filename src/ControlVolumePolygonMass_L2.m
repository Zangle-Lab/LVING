function [mask]=ControlVolumePolygonMass_L2(x, y, sz)


%% finding the points on sides of quadrilateral intersecting with the grid lines

for aa=1:length(x)

    xys{aa,1}=[x(aa),y(aa)];

    % for the x intercepts
    if aa<length(x)
        s=x(aa:aa+1);
        q=y(aa:aa+1);
    else
        s(1)=x(aa);
        s(2)=x(1);
        q(1)=y(aa);
        q(2)=y(1);
    end

    [x1,inl]=min(s);
    [x2,inh]=max(s);
    y1=q(inl);
    y2=q(inh);
    xi=ceil(x1);
    bb=1;
    cc=0;

    while (xi+cc)<x2  
        xys{aa,bb+1}=[(xi+cc),((((y2-y1)/(x2-x1))*(xi+cc-x1))+y1)];
        bb=bb+1;
        cc=cc+1;
    end
    
    % for the y intercepts
    if aa<length(x)
        s=y(aa:aa+1);
        q=x(aa:aa+1);
    else
        s(1)=y(aa);
        s(2)=y(1);
        q(1)=x(aa);
        q(2)=x(1);
    end
    [y1,inl]=min(s);
    [y2,inh]=max(s);
    x1=q(inl);
    x2=q(inh);
    cc=0;
    yi=ceil(y1);

    while (yi+cc)<y2  
        xys{aa,bb+1}=[((((x2-x1)/(y2-y1))*(yi+cc-y1))+x1),(yi+cc)];
        bb=bb+1;
        cc=cc+1;
    end

end

%% finding the grid points lying inside the quadrilateral
% upper and lower bound of quadrilateral on the grid

xl=min(x);
xh=max(x);
yl=min(y);
yh=max(y);
if ceil(yl)==yl
    ystart=yl+1;
else
    ystart=ceil(yl);
end
if floor(yh)==yh
    yend=yh-1;
else
    yend=floor(yh);
end
if ceil(xl)==xl
    xstart=xl+1;
else
    xstart=ceil(xl);
end
if floor(xh)==xh
    xend=xh-1;
else
    xend=floor(xh);
end
yint=ystart:yend;
xint=xstart:xend;

%test and record the points inside the quadrilateral

ss=size(xys);

for dd=1:length(xint)

    for ee=1:length(yint)
        Pointx=xint(dd);
        Pointy=yint(ee);
        pntnmb=0;
        Pnty=[];

        for ff=1:ss(1)
            for gg=1:ss(2)
                if isempty(xys{ff,gg})~=1
                    if xys{ff,gg}(1,1)==Pointx
                        if Pointy<=xys{ff,gg}(1,2)
                            pntnmb=pntnmb+1;
                            Pnty(pntnmb)=xys{ff,gg}(1,2);
                        end
                    end
                end
            end
        end

        QPnty=unique(Pnty);
        NumPnt=length(QPnty);

        if rem(NumPnt,2)~=0
            PntMtx{dd,ee}=[Pointx,Pointy];
        else
            PntMtx{dd,ee}=[];
        end

    end

end

%% finding the area inside each polygon formed on every element of grid
ygrid=floor(yl):yh;
xgrid=floor(xl):xh;

for hh=1:length(xgrid)

    for ii=1:length(ygrid)
        count=0;

        %list the points inside each pixel
        for jj=1:length(xint)
            for kk=1:length(yint)
                if isempty(PntMtx{jj,kk})~=1
                    if PntMtx{jj,kk}(1,1)<=(xgrid(hh)+1) && PntMtx{jj,kk}(1,1)>=xgrid(hh)
                        if PntMtx{jj,kk}(1,2)<=(ygrid(ii)+1) && PntMtx{jj,kk}(1,2)>=ygrid(ii)
                            count=count+1;
                            Vertice{hh,ii}(count,1:2)=PntMtx{jj,kk};
                        end
                    end
                end
            end
        end
        for jj=1:ss(1)
            for kk=1:ss(2)
                if isempty(xys{jj,kk})~=1
                    if xys{jj,kk}(1,1)<=(xgrid(hh)+1) && xys{jj,kk}(1,1)>=xgrid(hh)
                        if xys{jj,kk}(1,2)<=(ygrid(ii)+1) && xys{jj,kk}(1,2)>=ygrid(ii)
                            count=count+1;
                            Vertice{hh,ii}(count,1:2)=xys{jj,kk};
                        end
                    end
                end
            end
        end


        % input to shoelace formula requires listing vertices in clockwise order
        if count>0
            xv=Vertice{hh,ii}(:,1);
            yv=Vertice{hh,ii}(:,2);
            cx = mean(xv);
            cy = mean(yv);
            a = atan2(yv - cy, xv - cx);
            [~, order] = sort(a);
            OrderdedVertice{hh,ii}(:,1) = xv(order);
            OrderdedVertice{hh,ii}(:,2) = yv(order);
            FractionArea(hh,ii)=polyarea(OrderdedVertice{hh,ii}(:,1),OrderdedVertice{hh,ii}(:,2));
        end

    end

end

%% creating desired mask to be multiplied with the mass image
mask=zeros(sz(1),sz(2));
mask(floor(xl):floor(xl)+length(xgrid)-1,floor(yl):floor(yl)+length(ygrid)-1)=FractionArea;
mask=transpose(mask);