N = 200;
M = 5;
t=[0:N-1];
w=zeros(M,N); 
mu=0.2;%input('mu = ');
y(1) = 0.0;
y(2) = 0.0;
for j = 3:N
 y(j) = 0.95*y(j-1) - 0.195*y(j-2); 
end

x = y+randn(1,N)*0.5;
%x= y;
d = y;
for i=(M+1):N
   e(i) = d(i) -  x((i-(M)+1):i)*w(:,i);
   w(:,i+1) = w(:,i) + mu * e(i) * x((i-(M)+1):i)';
end
for i=(M+1):N
    yd(i) = x((i-(M)+1):i)*w(:,i);  
end
  subplot(221),plot(t,d),ylabel('Desired Signal'),
  subplot(222),plot(t,x),ylabel('Input Signal+Noise'),
  subplot(223),plot(t,e),ylabel('Error'),
  subplot(224),plot(t,yd),ylabel('Adaptive Desired output')