
% function to check if the polygon is convex

function [check]=checkConvex(x,y)

ds1(1)=x(2)-x(1);
ds1(2)=x(3)-x(2);
ds1(3)=y(2)-y(1);
ds1(4)=y(3)-y(2);
CrossPdt(1)=(ds1(1)*ds1(4))-(ds1(2)*ds1(3));

ds2(1)=x(3)-x(2);
ds2(2)=x(4)-x(3);
ds2(3)=y(3)-y(2);
ds2(4)=y(4)-y(3);
CrossPdt(2)=(ds2(1)*ds2(4))-(ds2(2)*ds2(3));

ds3(1)=x(4)-x(3);
ds3(2)=x(1)-x(4);
ds3(3)=y(4)-y(3);
ds3(4)=y(1)-y(4);
CrossPdt(3)=(ds3(1)*ds3(4))-(ds3(2)*ds3(3));

ds4(1)=x(1)-x(4);
ds4(2)=x(2)-x(1);
ds4(3)=y(1)-y(4);
ds4(4)=y(2)-y(1);
CrossPdt(4)=(ds4(1)*ds4(4))-(ds4(2)*ds4(3));

E=CrossPdt>0;
F=CrossPdt<0;

if sum(E)==4 || sum(F)==4
    check=1;
else
    check=0;
end