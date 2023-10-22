
%function to find the peak of A, using subpixel grid Xsub,Ysub, of size
% uses fast gaussian ps function from matlab file exchange to find peak

function [ip, jp,results] = findvalley_v3(A)

sz = size(A);
[row,col] = find(ismember((A), min(min(A(:)))));
if length(row)==1 && length(col)==1
    if row>1 && row<sz(1)-1
        if col>1 && col<sz(2)-1
            E=A(row-1:row+1,col-1:col+1);
            a=max(E(:));
            Q=double(a-E);
            results = psfFit(Q);   % find inside 'private' folder
            ip=(-((sz(1)+1)/2)+col)-(2-results(1));
            jp=(-((sz(2)+1)/2)+row)-(2-results(2));
        else
            ip=0;
            jp=0;
        end
    else
        ip=0;
        jp=0;
    end
else
    ip=0;
    jp=0;
end