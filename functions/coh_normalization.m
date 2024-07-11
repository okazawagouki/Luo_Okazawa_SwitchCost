function coh = coh_normalization(coh, cond)
%   coherence normalization

%   dots task:
%   0   0.032   0.064   0.128   0.256   0.512
%==>0   0.064   0.128   0.256   0.512   1      so as to match the color task

%   color task:
%   0           0.064   0.128   0.256   0.512   1

coh(cond==1) = coh(cond==1)*2; % cond1 is dots task
coh(cond==1 & coh==-1.024) = -1;
coh(cond==1 & coh==1.024) = 1;
end


















