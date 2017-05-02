function [ code ] = Code_Generator( system, sample_1_2, svnum, N)
% if nargin 
if(sample_1_2 == 1)
    sample = sample_1_2;
else
    sample = 1;
end

switch system
    case {'gps','sbas'}
        % Generate Gold_Code
        g2s = [5; 6; 7; 8; 17; 18; 139; 140; 141; 251; 252; 254; 255; 256; 257; 258;...
               469; 470; 471; 472; 473; 474; 509; 512; 513; 514; 515; 516; 859; 860;...
               861; 862;... %% 1- 32 for GPS
               145; 175; 52; 21; 237; 235; 886; 657; 634; 762; 355; 1012; 176; 603;... 
               130; 359; 595; 68; 386]; %% 120 - 138 for SBAS

        if sample < 1
            sample = 1;
        end

        code = zeros(1023*sample, 1);

        g1      = zeros(1, 1023);
        g2      = zeros(1, 1023);

        if (svnum < 0) || (svnum > 32 && svnum < 120) || (svnum > 138)
            return;
        end

        if svnum < 33
            codeNum = svnum;
        else
            codeNum = svnum - 120 + 33;
        end

        g2shift=g2s(codeNum, 1);
        % ******* Generate G1 code *******
        % load shift register
        reg = -1*ones(1,10);
        for i = 1:1023,
            g1(i) = reg(10);
            save1 = reg(3)*reg(10);
            reg(1,2:10) = reg(1:1:9);
            reg(1) = save1;
        end,
        % ******* Generate G2 code *******
        % load shift register
        reg = -1*ones(1,10);
        for i = 1:1023,
            g2(i) = reg(10);
            save2 = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
            reg(1,2:10) = reg(1:1:9);
            reg(1) = save2;
        end,
        % ******* Shift G2 code *******
        g2tmp(1,1:g2shift)=g2(1,1023-g2shift+1:1023);
        g2tmp(1,g2shift+1:1023)=g2(1,1:1023-g2shift);
        g2 = g2tmp;
        % ******* Form single sample C/A code by multiplying G1 and G2
        ss_ca = g1.*g2;
        code_used=-ss_ca;

        for j = 0:sample-1
         code(1+j:sample:end, 1) = code_used(:);
        end
        
    case 'gps_L5_I'
        [isI, isQ] = InitialState_L5_B(svnum);
        
        code_used = L5_I5_gen(svnum, 1, isI);
        code = zeros(10230*sample, 1);
        for j = 0:sample-1
            code(1+j:sample:end, 1) = code_used(:);
        end
        
    case 'gps_L5_Q'
        [isI, isQ] = InitialState_L5_B(svnum);
        
        code_used = L5_Q5_gen(svnum, 1, isQ);
        
        code = zeros(10230*sample, 1);
        for j = 0:sample-1
            code(1+j:sample:end, 1) = code_used(:);
        end
        
        
    case 'glonass'
        %Generate M_code
        if nargin == 1
             sample = 1;
        end

        if sample < 1
            sample = 1;
        end

        code_used = zeros(511, 1);
        code      = zeros(511*sample, 1);

        shift_reg=ones(9, 1);

        for i=1:511 
            code_used(i, 1) = shift_reg(7, 1); 
            mod2 = xor(shift_reg(5, 1), shift_reg(9, 1));
            shift_reg(2:9, 1) = shift_reg(1:8, 1);
            shift_reg(1, 1) = mod2;
        end
        
        

        for j = 0:sample-1
         code(1+j:sample:end, 1) = -2*code_used(:)+1;
        end
        
    case 'E1_B_primary'
        load E1_B_bin.mat;
        
        if svnum < 51
            codeNum = svnum;
        else
            codeNum = 1;
        end

        
        code = zeros(4092*sample, 1);
        code_used = (-1).^ E1_B_bin(codeNum, :);
        
        for j = 0:sample-1
            code((1+j):sample:end, 1) = code_used(:);
        end
        
    case 'E1_C_primary'
        load E1_C_bin.mat;
        
        if svnum < 51
            codeNum = svnum;
        else
            codeNum = 1;
        end

        
        code = zeros(4092*sample, 1);
        code_used = (-1).^ E1_C_bin(codeNum, :);
        
        for j = 0:sample-1
            code((1+j):sample:end, 1) = code_used(:);
        end
        
    case 'E1_C_secondary'
        load secondary_bin.mat;
        if svnum < 50
            codeNum = 50+svnum;
        else
            codeNum = 51;
        end

        code = (-1).^secondary_codes_bin.CS100{codeNum, :};
        return

    case 'bds' %  »“ј…— јя ЅЁ…ƒ”ќ она же  ќћѕј——
        
        svShift = [
            1 3 ;%//1
            1 4 ;%//1
            1 5 ;%//2
            1 6 ;%//3
            1 8 ;%//4
            1 9 ;%//5
            1 10;%//6
            1 11;%//7
            2 7 ;%//8
            3 4 ;%//9
            3 5 ;%//10
            3 6 ;%//11
            3 8 ;%//12
            3 9 ;%//13
            3 10;%//14
            3 11;%//15
            4 5 ;%//16
            4 6 ;%//17
            4 8 ;%//18
            4 9 ;%//19
            4 10;%//20
            4 11;%//21
            5 6 ;%//22
            5 8 ;%//23
            5 9 ;%//24
            5 10;%//25
            5 11;%//26
            6 8;%//27
            6 9 ;%//28
            6 10;%//29
            6 11;%//30
            8 9 ;%//31
            8 10;%//32
            8 11;%//33
            9 10;%//34
            9 11;%//35
            10 11;%//37
        ];
    
        %Generate Gold_Code
        if nargin == 1
             sample = 1;
        end

        if sample < 1
            sample = 1;
        end
        
        codeLen = 2046;
        code = zeros(codeLen*sample, 1);
        
        for increment = 1:sample
            for inc = 1:length(svnum)
                sv_1 = svShift(svnum(inc), 1);
                sv_2 = svShift(svnum(inc), 2); 
                reg_G1 = ones(11, 1);
                reg_G2 = ones(11, 1);
                for i = 1:11
                    reg_G1(i, 1) = mod(i+1, 2);
                    reg_G2(i, 1) = mod(i+1, 2);
                end
                for i = 1:codeLen
                    if ( mod( reg_G1(11) + reg_G2(sv_1) + reg_G2(sv_2), 2) == 1)
                        code(i + (increment - 1)*1023) = - 1;
                    else
                        code(i + (increment - 1)*1023) = 1;
                    end;
                    temp_1 = mod( reg_G1(1) + reg_G1(7) + reg_G1(8) + reg_G1(9) + reg_G1(10) + reg_G1(11) , 2);
                    temp_2 = mod( reg_G2(1) + reg_G2(2) + reg_G2(3) + reg_G2(4)...
                                + reg_G2(5) + reg_G2(8) + reg_G2(9) + reg_G2(11), 2);
                    reg_G1 = circshift(reg_G1, 1);
                    reg_G2 = circshift(reg_G2, 1);
                    reg_G1(1) = temp_1;
                    reg_G2(1) = temp_2;
                end
            end
        end
        
    case 'random_N'
        code = zeros(N*sample, 1);
        code0 = (-1).^(random('unid', 100, N, 1));
        for increment = 0:sample-1
            code((1+increment):2:N*sample) = code0(:);
        end
        
    case 'random_511'
        code = zeros(511*sample, 1);
        code0 = (-1).^(random('unid', 100, 511, 1));
        for increment = 0:sample-1
            code((1+increment):2:511*sample) = code0(:);
        end
        
    case 'random_1023'
        code = zeros(1023*sample, 1);
        code0 = (-1).^(random('unid', 100, 1023, 1));
        for increment = 0:sample-1
            code((1+increment):2:1023*sample) = code0(:);
        end
        
    case 'random_10230'
        code = zeros(10230*sample, 1);
        code0 = (-1).^(random('unid', 100, 10230, 1));
        for increment = 0:sample-1
            code((1+increment):2:10230*sample) = code0(:);
        end
        
    case 'glonass_all_band'
        %Generate M_code
        if nargin == 1
             sample = 1;
        end

        if sample < 1
            sample = 1;
        end

        code_used = zeros(511, 1);
        code      = zeros(511*sample, 1);

        shift_reg=ones(9, 1);

        for i=1:511 
            code_used(i, 1) = shift_reg(7, 1); 
            mod2 = xor(shift_reg(5, 1), shift_reg(9, 1));
            shift_reg(2:9, 1) = shift_reg(1:8, 1);
            shift_reg(1, 1) = mod2;
        end

        for j = 0:sample-1
            code(1+j:sample:end, 1) = -2*code_used(:)+1;
        end
        
    otherwise
        disp('—истема с таким именем не используетс€ в данном проекте!');
        code = zeros(sample_1_2(1), 1);
end

if(sample_1_2 > 1)
    freaq2 = sample_1_2;
    freaq1 = sample_1_2;
    code = RawResampler( code, freaq1, freaq2);
end

if(strcmp(system, 'glonass_all_band'))
    % фильтруем
    Hd = Filter_Glonass_code_40_clc;
    code = filter(Hd, code);
    
%     сдвигаем на номер литеры
%     prnFrq = 511e3;
%     deltaT = sample_1_2(2)/(prnFrq*sample_1_2(1));
%     delta_F = 562500; % Hz
%     freq_Shift = delta_F*(svnum-8);
%     freq_Shift
%     code = Mixer(code, freq_Shift, deltaT);
    
%     figure(123)
%     spectr = abs(fftshift(fft(code)));
%     plot(20*log10(spectr/max(spectr)));
%     sum(spectr)
%     hold all
end

end

