f = NaN(30,1);
f(1) = mean_ft(4);
c_ = OutParams(8:13);
A_ = OutParams(14:19);
B_ = OutParams(20:25);

for i=2:30
    f(i) = c_(4)+A_(4)*f(i-1)+B_(4)*Res.score(4,i-1);
end