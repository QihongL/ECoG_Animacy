A=[3 4 5 6 7];
B=[4 6 7];
[tf,loc]=ismember(A,B);
idx=[1:length(A)];
idx=idx(tf);
idx=idx(loc(tf));
disp(A(idx))