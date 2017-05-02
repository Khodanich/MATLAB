function  Q  = Quality(Corr_massive,lengthCorr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MaxCorr = max(Corr_massive);

for i = 1:lengthCorr
    if  Corr_massive(i) == MaxCorr
        MaxNum = i;
         break;
    end
end
MaxNumL = MaxNum - 1;
MaxNumR = MaxNum + 1;
if MaxNum == 1
    MeanCorr = mean(Corr_massive) - (Corr_massive(MaxNum) + Corr_massive(MaxNumR))/lengthCorr;
elseif MaxNum == 4092    
    MeanCorr = mean(Corr_massive) - (Corr_massive(MaxNum) + Corr_massive(MaxNumL))/lengthCorr;
else 
    MeanCorr = mean(Corr_massive) - (Corr_massive(MaxNum) + Corr_massive(MaxNumL) + Corr_massive(MaxNumR))/lengthCorr;
end

Q = (MaxCorr - abs(MeanCorr))/MaxCorr * 100;

end

