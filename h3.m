mu = 0.2;
M = 5;  %порядок фильтра
I = 1000;%
Input = zeros(1,I);
Desired = zeros(1,I);

W = zeros(1,M);
X = zeros(1,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [0.0625, 0.125, 0.5, 0.25, 1];

for i = 1:I
    Input = rand(1 , I);    
end
for i = 1:I
    for j = 1:M
        if ((i - j) > 0)
            Desired(i) = Desired(i) + Input(i - j)*H(j);
        end
    end
end

for T = 1:I
    
    for m = T:-1:(T-M)  % подумать над диапазоном значений
        if (m >= 0)
           X(M + (m - T)) = Input(m);         
        else break;
        end
    end
    
    D = Desired(T); % желаемый сигнал
    Y = 0;
    for i = 1:M
       Y = Y + W(i)*X(i);   % 
    end   
    
    E = D - Y;  %ошибка
    
    for i = 1:M
        W(i) = W(i) + mu*E*X(i); %обновление весового коэффициента
    end    
end
plot(1:I,E);