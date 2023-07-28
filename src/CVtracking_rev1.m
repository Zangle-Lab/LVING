function [SOPL2]=CVtracking_rev1(xex,yey,xcg,SOPL)
sz=size(SOPL);
sez=sz(1);
SOPL1=zeros(sz(1),sz(2),'single');
parfor pp=xcg+1:sez-xcg+1
    for qq=xcg+1:sez-xcg+1
        xs1=[xex(pp-(xcg/2),qq-(xcg/2)),xex(pp+(xcg/2),qq-(xcg/2)),xex(pp+(xcg/2),qq+(xcg/2)),xex(pp-(xcg/2),qq+(xcg/2))];
        ys1=[yey(pp-(xcg/2),qq-(xcg/2)),yey(pp+(xcg/2),qq-(xcg/2)),yey(pp+(xcg/2),qq+(xcg/2)),yey(pp-(xcg/2),qq+(xcg/2))];
        SOPL1(pp,qq)=ControlVolumePolygonMass_rev6(xs1,ys1,SOPL);
%         SOPL1(pp,qq)=sum(SOPL.*mask);
    end
end
SOPL2=(SOPL1(xcg+1:sz(1)-xcg,xcg+1:sz(2)-xcg))/(xcg*xcg);
