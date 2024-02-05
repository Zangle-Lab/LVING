
% This function determines if the tracking of control volume has resulted
% in a convex or concave polygon. The chnage in shape of the control volume
% decides the math for computing the area of the polygon

function [Value]=ControlVolumePolygonMass(x, y,I)

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

    EF=CrossPdt<0; 
    %if sum of EF is 0:convex.. 1:concave-not intersecting.. 2:concave-intersecting
    
    % if polygon intersects on itself, find the intersection point and confirm
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

        % if 1-2 and 3-4 are intersecting
        if DotPdt>0 && DotPdt<LngSq 
            xx1=[x(1),X1,x(4)];
            yy1=[y(1),Y1,y(4)];
            Clockdir=(yy1(2)-yy1(1))*(xx1(3)-xx1(2))-(yy1(3)-yy1(2))*(xx1(2)-xx1(1));
            [mask1]=ControlVolumePolygonMass_L2(xx1,yy1,sz);
            Value1=sum(sum(I.*mask1));
            xx2=[X1,x(2),x(3)];
            yy2=[Y1,y(2),y(3)];
            [mask2]=ControlVolumePolygonMass_L2(xx2,yy2,sz);
            Value2=sum(sum(I.*mask2));
            if Clockdir>0
                Value=Value1-Value2;
            else
                Value=Value2-Value1;
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

        % if 2-3 and 4-1 are intersecting
        if DotPdt>0 && DotPdt<LngSq  
            xx1=[x(1),x(2),X2];
            yy1=[y(1),y(2),Y2];
            Clockdir=(yy1(2)-yy1(1))*(xx1(3)-xx1(2))-(yy1(3)-yy1(2))*(xx1(2)-xx1(1));
            [mask1]=ControlVolumePolygonMass_L2(xx1,yy1,sz);
            Value1=sum(sum(I.*mask1));
            xx2=[X2,x(3),x(4)];
            yy2=[Y2,y(3),y(4)];
            [mask2]=ControlVolumePolygonMass_L2(xx2,yy2,sz);
            Value2=sum(sum(I.*mask2));
            if Clockdir>0
                Value=Value1-Value2;
            else
                Value=Value2-Value1; 
            end
        end
    
    else

        Clockdir=(y(2)-y(1))*(x(3)-x(2))-(y(3)-y(2))*(x(2)-x(1));
        if Clockdir>0
            [mask]=ControlVolumePolygonMass_L2(x,y,sz);
            Value=sum(sum(I.*mask));
        else
            [maskb]=ControlVolumePolygonMass_L2(x,y,sz);
            Valueb=sum(sum(I.*maskb));
            Value=-Valueb; mask=-maskb;
        end
    
    end

else
    
    Value=0;

end