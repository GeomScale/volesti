function [C,b,c0]=getEllMat(array,d,nrows)

C=zeros(d,d);
b=ones(d,1)/(d+1);
c0=0.00002;
EllMat=zeros(nrows,d);
x_med=zeros(d,1);
for i=1:nrows
    for j=1:d
        EllMat(i,j)=array( (i-1)*d+j );
    end
end

for i=1:nrows
    x=EllMat(i,:);
    C=C+x'*x;
    x_med=x_med+x';
end

C=C/nrows; x_med=x_med/nrows;

C=C-x_med*x_med';
end