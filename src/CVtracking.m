
% function to compute the deform phase image relative to prior 
% phase image using calculated grid deformation
% SOPL is current phase and SOPL2 is after deformation
% xex and yey defines grid deformation defined on grid of mesh size xcg

function [SOPL2]=CVtracking(xex,yey,xcg,SOPL)
sz=size(SOPL);
sez=sz(1);
SOPL1=zeros(sz(1),sz(2),'single');
parfor pp=xcg+1:sez-xcg+1
    for qq=xcg+1:sez-xcg+1
        xs1=[xex(pp-(xcg/2),qq-(xcg/2)),xex(pp+(xcg/2),qq-(xcg/2)),xex(pp+(xcg/2),qq+(xcg/2)),xex(pp-(xcg/2),qq+(xcg/2))];
        ys1=[yey(pp-(xcg/2),qq-(xcg/2)),yey(pp+(xcg/2),qq-(xcg/2)),yey(pp+(xcg/2),qq+(xcg/2)),yey(pp-(xcg/2),qq+(xcg/2))];
        SOPL1(pp,qq)=ControlVolumePolygonMass(xs1,ys1,SOPL);
    end
end
SOPL2=(SOPL1(xcg+1:sz(1)-xcg,xcg+1:sz(2)-xcg))/(xcg*xcg);
