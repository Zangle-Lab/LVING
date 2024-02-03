%% LVING : figures showing the proof of relevance of the assumptions used in LVING computation

% path to the dataset input
load('K:\Data\Soorya\RPEDrugTest_2020\Pn7_Ethanol0.04uL_21May2020\MassGenResults_rev136\Results_2hr\WS1_cell16.mat')

% compute the elements used in the CV computation
for pp=8:492
    for qq=8:492

        delMt = sum(sum((Abkg_mass(pp:pp+4,qq:qq+4,2)-Abkg_mass(pp:pp+4,qq:qq+4,1))));  % vel in pixel/min

        delX = (dX(pp,qq,1).*Abkg_mass(pp,qq,1)) - (dX(pp,qq+4,1).*Abkg_mass(pp,qq+4,1)) ;

        delY = (dY(pp,qq,1).*Abkg_mass(pp,qq,1)) - (dY(pp+4,qq,1).*Abkg_mass(pp+4,qq,1)) ;

        delR(pp,qq) = delMt-(delX+delY);
        
        delM(pp,qq)=delMt;
        
        delD(pp,qq) = delX+delY ;
        
    end
end

%% bar plot of DM/Dt, DM/Dt + vx.delx + vy.dely , vx.delx + vy.dely
figure(1); 
h=bar([mean(delM(delM~=0)),mean(delR(delR~=0)), mean(delD(delD~=0))]);
set(gca, 'XTickLabel', {'growth','growth+grad(vm)','grad(vm)'});
hold on;
errorbar([1,2,3],[mean(delM(delM~=0)),mean(delR(delR~=0)), mean(delD(delD~=0))],...
    [std(delM(delM~=0))/sqrt(nnz(delM)),std(delR(delR~=0))/sqrt(nnz(delR)),std(delD(delD~=0))/sqrt(nnz(delD))],...
    [std(delM(delM~=0))/sqrt(nnz(delM)),std(delR(delR~=0))/sqrt(nnz(delR)),std(delD(delD~=0))/sqrt(nnz(delD))],...
    '.','Color','b');
hold on;
% ydata=[nonzeros(delM),nonzeros(delR)]; [r, c] = size(ydata); xdata = repmat(1:c, r, 1);
% scatter(xdata(:), ydata(:), 'k.', 'jitter','on', 'jitterAmount', 0.05);
set(gcf,'Color','w'); ylim([-0.00005 0.0002]);
ylabel('CV change in mass (pg/min)');

%% histogram of grad(vm) overlapped with DM/Dt
figure(2); 
histogram(delD);
hold on; histogram(delM);
xlim([-0.005 0.005]); set(gcf,'Color','w');


