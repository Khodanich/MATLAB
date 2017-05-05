mu = 0.01;
%M = 5;  %порядок фильтра
M = 5;
I = 10000;%
Input = zeros(1,I);
Desired = zeros(1,I);
%NoiseDb = 100;
w0=0.001;  phi=0.1;

W = zeros(1,M);
X = zeros(1,M); 


MainFilt = [1, 0.5, 0.25, 0.125, 0.0625];
H = fliplr(MainFilt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H = 1:5;

% for i = 1:I
%     Input = rand(1 ,I);  
% end
Input1 = sin(2*pi*[1:I]*w0+phi);
Input = sin(2*pi*[1:I]*w0+phi) + 2*rand(1,I);

for i = 1:I
    for j = 1:M
        if ((i - j) > 0)
            Desired(i) = Desired(i) + Input(i - j)*H(j);
        end
    end
end

for T = 1:I
    
    for m = T:-1:(T - M)
        if (m >= 1)
           X(M + (m - T) + 1) = Input(m);         
        else break;
        end
    end
    
    D = Desired(T); % желаемый сигнал
    Y = 0;
    for i = 1:M
       Y = Y + W(i)*X(i);   % 
    end   
    
    E = D - Y;  %ошибка
    Error(T) = E;
    for i = 1:M
        W(i) = W(i) + mu*E*X(i); %обновление весового коэффициента
    end    
end
plot(1:I,Error);



for i = 1:I
    for j = 1:M
        if ((i - j) > 0)
            Desired(i) = Desired(i) + Input(i - j)*X(j);
        end
    end
end
%plot(1:I,Desired);










