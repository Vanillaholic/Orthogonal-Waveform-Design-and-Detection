function [A, vel] = wideband_ambiguity(x1, x2, r, b, len, c, v_range,r_range, fs)


    rs_bsig1 = resample(x1,1,1);        %对信号不进行重采样
    % delay
    N = size(rs_bsig1,2);               %获取重采样的信号大小
    delay = round((r_range/c)*fs);         %TO DO:获取1000m距离
    obsvN = delay+N;                    %obsvN 总的观测的长度
    s_n = [rs_bsig1 zeros(1,obsvN-N)];  %s_n 在原始信号后面补零
    sig = s_n(1:obsvN-delay);           %截取的信号
    bsig_no = [zeros(1,delay) sig];     %对sig添加延迟

    clear obsvN sig N;

    vel_del = c/(r*(b^len-1));          % vel_del 速度分辨率
    vel = 0:vel_del:v_range;               
    vel = [-vel(end:-1:2) vel];         % vel 速度范围：其中存储的是不同的速度，包括正负值
    eta = 1+(vel/c);                    %获取(1+v/c)的分子和分母
    [p,q]= rat(eta);                    

    ambig1 = cell(1, length(vel));      %初始化两个cell数组


    for i =1:length(vel)  % 遍历vel速度数组
        re_samp_bsig1 = resample(x2,p(i),q(i));     % p>q expansion and p<q compression
    
    % --------在时域上进行相关处理----------------------------------- 
        if length(re_samp_bsig1)>length(bsig_no)   %如果重采样后的信号长度大于延迟信号，补零以匹配长度。
            na = length(re_samp_bsig1)-length(bsig_no);
            bsig_no= [bsig_no zeros(1,na+1)];
        end
        ambig1{i}=  abs(matchFilter(bsig_no,re_samp_bsig1,'none'));
        nl = length(bsig_no);
        ambigh1(1:(nl),i) =(ambig1{i});
    end

    Max1= max(max(ambigh1));
    A=(abs(ambigh1 ./ Max1));

%重采样模糊函数
%% resample the ambiguity funnction;
    for i =1:length(vel)
        re_ambigh1(:,i)=resample(A(:,i),1,6.);

    end
    A=re_ambigh1;
end


