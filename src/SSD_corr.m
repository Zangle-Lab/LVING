
% function returns the SSD computed between two images (CurrD and NextD) 
% over overlapping regions of size grid
% SS is ndarray of all SSDs computed 

function [SS]=SSD_corr(CurrD,NextD,grid)

[l,c]=size(CurrD);
X=ones(l+(2*grid),c+(2*grid))*mean2(CurrD);
X(grid+1:l+grid,grid+1:c+grid)=CurrD;

SS=zeros((2*grid)+1,(2*grid)+1,l,c,'single','gpuArray');
for aa=1:(2*grid)+1
    for bb=1:(2*grid)+1
        SS(aa,bb,:,:)=movsum(movsum(((X(aa:l+aa-1,bb:c+bb-1)-NextD).^2),grid,1),grid,2);
    end
end