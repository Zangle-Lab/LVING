function [Value]=ControlVolumePolygonMass_rev6(x, y,I)

%% size of image to match with the size of mask, all x are row no and y is column no
sz=size(I);
swx=x>=sz(1)-2;
bwx=x<=1;
sbx=sum(swx(:))+sum(bwx(:));
swy=y>=sz(2)-2;
bwy=y<=1;
sby=sum(swy(:))+sum(bwy(:));
if sbx==0 && sby==0

    %% layout details of the polygon
    % is the polygon concave or twisted
    [CrossPdt]=checkConvex(x,y);
    EF=CrossPdt<0; %if sum of EF is 0:convex.. 1:concave-not intersecting.. 2:concave-intersecting
    % if polygon intersects on itself, find the intersection point and confirm
%     InterPnt=[];
    if sum(EF)==2
        % if 1-2 and 3-4 are intersecting
        if x(1)~=x(2) && x(3)~=x(4)
            X1=(y(1)-y(3)+(x(3)*((y(4)-y(3))/(x(4)-x(3))))-(x(1)*((y(2)-y(1))/(x(2)-x(1)))))/(((y(4)-y(3))/(x(4)-x(3)))-((y(2)-y(1))/(x(2)-x(1))));
            Y1=(((y(2)-y(1))/(x(2)-x(1)))*X1)+y(1)-(x(1)*((y(2)-y(1))/(x(2)-x(1))));
        end
        if x(1)==x(2)
            X1=x(1);
            Y1=((X1-x(3))*((y(4)-y(3))/(x(4)-x(3))))+y(3);
        end
        if x(3)==x(4)
            X1=x(3);
            Y1=((X1-x(1))*((y(2)-y(1))/(x(2)-x(1))))+y(1);
        end
        DotPdt=(X1 -x(1)) * (x(2) - x(1)) + (Y1 - y(1))*(y(2) - y(1));
        LngSq=(x(2) - x(1))*(x(2) - x(1)) + (y(2) - y(1))*(y(2) - y(1));
        if DotPdt>0 && DotPdt<LngSq % if 1-2 and 3-4 are intersecting
%             InterPnt=[X1,Y1];
            xx1=[x(1),X1,x(4)];
            yy1=[y(1),Y1,y(4)];
            Clockdir=(yy1(2)-yy1(1))*(xx1(3)-xx1(2))-(yy1(3)-yy1(2))*(xx1(2)-xx1(1));
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(xx1,yy1,sz);
            [mask1]=ControlVolumePolygonMass_revL2(xx1,yy1,sz);
            Value1=sum(sum(I.*mask1));
%             Value1=sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
            xx2=[X1,x(2),x(3)];
            yy2=[Y1,y(2),y(3)];
            [mask2]=ControlVolumePolygonMass_revL2(xx2,yy2,sz);
            Value2=sum(sum(I.*mask2));
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(xx2,yy2,sz);
%             Value2=sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
            if Clockdir>0
                Value=Value1-Value2;
%                 mask=mask1-mask2;
            else
                Value=Value2-Value1;
%                 mask=mask2-mask1;
            end
        end

        % if 2-3 and 4-1 are intersecting
        if x(2)~=x(3) && x(4)~=x(1) 
            X2=(y(2)-y(4)+(x(4)*((y(1)-y(4))/(x(1)-x(4))))-(x(2)*((y(3)-y(2))/(x(3)-x(2)))))/(((y(1)-y(4))/(x(1)-x(4)))-((y(3)-y(2))/(x(3)-x(2))));
            Y2=(((y(3)-y(2))/(x(3)-x(2)))*X2)+y(2)-(x(2)*((y(3)-y(2))/(x(3)-x(2))));
        end
        if x(2)==x(3)
            X2=x(2);
            Y2=((X2-x(4))*((y(1)-y(4))/(x(1)-x(4))))+y(4);
        end
        if x(4)==x(1)
            X2=x(4);
            Y2=((X2-x(2))*((y(3)-y(2))/(x(3)-x(2))))+y(2);
        end
        DotPdt=(X2 -x(2)) * (x(3) - x(2)) + (Y2 - y(2))*(y(3) - y(2));
        LngSq=(x(3) - x(2))*(x(3) - x(2)) + (y(3) - y(2))*(y(3) - y(2));
        if DotPdt>0 && DotPdt<LngSq  % if 2-3 and 4-1 are intersecting
%             InterPnt=[X2,Y2];
            xx1=[x(1),x(2),X2];
            yy1=[y(1),y(2),Y2];
            Clockdir=(yy1(2)-yy1(1))*(xx1(3)-xx1(2))-(yy1(3)-yy1(2))*(xx1(2)-xx1(1));
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(xx1,yy1,sz);
%             Value1=sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
            [mask1]=ControlVolumePolygonMass_revL2(xx1,yy1,sz);
            Value1=sum(sum(I.*mask1));
            xx2=[X2,x(3),x(4)];
            yy2=[Y2,y(3),y(4)];
            [mask2]=ControlVolumePolygonMass_revL2(xx2,yy2,sz);
            Value2=sum(sum(I.*mask2));
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(xx2,yy2,sz);
%             Value2=sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
            if Clockdir>0
                Value=Value1-Value2;
%                 mask=mask1-mask2;
            else
                Value=Value2-Value1; 
%                 mask=mask2-mask1;
            end
        end
    else
        Clockdir=(y(2)-y(1))*(x(3)-x(2))-(y(3)-y(2))*(x(2)-x(1));
        if Clockdir>0
            [mask]=ControlVolumePolygonMass_revL2(x,y,sz);
            Value=sum(sum(I.*mask));
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(x,y,sz);
%             Value=sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
        else
            [maskb]=ControlVolumePolygonMass_revL2(x,y,sz);
            Valueb=sum(sum(I.*maskb));
            Value=-Valueb; mask=-maskb;
%             [CoPnt,FractionArea]=ControlVolumePolygonMass_revL3(x,y,sz);
%             Value=-sum(sum(I(CoPnt(1):CoPnt(2),CoPnt(3):CoPnt(4)).*transpose(FractionArea)));
%             Value=-Valueb; 
%             mask=-maskb;
        end
    end
else
    Value=0;
%     mask=zeros(sz(1),sz(2),'single');
end

% totalArea=polygonArea(y,x,4);