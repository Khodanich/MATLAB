N=input('length of sequence N = ');
    t=[0:N-1];
   w0=0.001;  phi=0.1;
   d=sin(2*pi*[1:N]*w0+phi);
   x=d+randn(1,N)*0.5;
    w=zeros(1,N); 
   mu=input('mu = ');
   for i=1:N
   yd(i)=w*x'; 
  e(i) = d(i) - w * x';
 for m=1:N
 w(m) = w(m) + mu * e(i) * x(m);
  end
 end
    subplot(221),plot(t,d),ylabel('Desired Signal'),
  sub plot(222),plot(t,x),ylabel('Input Signal+Noise'),
  subplot(223),plot(t,e),ylabel('Error'),
  subplot(224),plot(t,yd),ylabel('Adaptive Desired output')