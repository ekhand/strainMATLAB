function [ x,y ] = interpolate_cscvn( f )
s=0;
x=zeros(3000,1);
y=zeros(3000,1);

for i=1:f.pieces 
    for t=0:0.5:abs(f.breaks(1,i)-f.breaks(1,i+1))
        s=s+1;
        x(s,:)=f.coefs(2*i-1,1)*t^3+f.coefs(2*i-1,2)*t^2+f.coefs(2*i-1,3)*t^1+f.coefs(2*i-1,4)*t^0;
        y(s,:)=f.coefs(2*i,1)*t^3+f.coefs(2*i,2)*t^2+f.coefs(2*i,3)*t^1+f.coefs(2*i,4)*t^0;
    end
    
    if i==f.pieces && t~=abs(f.breaks(1,i)-f.breaks(1,i+1))
        t=abs(f.breaks(1,i)-f.breaks(1,i+1));
        s=s+1;
        x(s,:)=f.coefs(2*i-1,1)*t^3+f.coefs(2*i-1,2)*t^2+f.coefs(2*i-1,3)*t^1+f.coefs(2*i-1,4)*t^0;
        y(s,:)=f.coefs(2*i,1)*t^3+f.coefs(2*i,2)*t^2+f.coefs(2*i,3)*t^1+f.coefs(2*i,4)*t^0;
    end
end
x=x(1:s,1);
y=y(1:s,1); 
end

% https://www.mathworks.com/matlabcentral/answers/406920-points-of-cscvn-function