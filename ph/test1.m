x=-1:0.1:1;              % ���������� ����� �� ��� Ox 
y=-2:0.1:2;              % ���������� ����� �� ��� Oy

[X,Y]=meshgrid(x,y);
Z=exp(-X.^2-Y.^2); 
plot3(X,Y,Z);
surf(X,Y,Z); 