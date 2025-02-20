% function [AF, tauAxis, vAxis] = widebandAF(x, fs, c,vAxis, maxDelay)
% % WIDEBANDAF 计算宽带模糊函数 (Wideband Ambiguity Function)
% %
% %   [AF, tauAxis, vAxis] = widebandAF(x, fs, vAxis, c, maxDelay)
% %
% % 输入：
% %   x       : 输入信号 (行向量或列向量均可)
% %   fs      : 采样率 (Hz)
% %   vAxis   : 速度向量，例如 linspace(-20, 20, 41)
% %   c       : 波速，例如声学场景中 c = 1500 m/s，电磁波中 c = 3e8 m/s
% %   maxDelay: 最大延迟范围（以采样点计），默认为 length(x)-1
% %
% % 输出：
% %   AF      : 模糊函数幅度值 (二维矩阵, 行对应延迟, 列对应速度)
% %   tauAxis : 延迟对应的时间轴 (秒)
% %   vAxis   : 速度轴 (与输入相同)
% %
% % 原理：
% %   对于每个速度 v，计算多普勒伸缩因子 alpha = 1 + v/c，
% %   对信号进行重采样来近似多普勒效应，然后与原始信号在各个延迟下进行匹配滤波或相关，
% %   结果存储在 AF 矩阵中。
% 
%     % 保证输入信号为行向量
%     x = x(:).';  
%     Nx = length(x);
% 
%     % 如果没有指定 maxDelay，则默认取全信号范围
%     if nargin < 5 || isempty(maxDelay)
%         maxDelay = Nx - 1;
%     end
% 
%     % 定义延迟范围，范围为 [-maxDelay, maxDelay]
%     tauSamples = -maxDelay:maxDelay;
%     Ntau = length(tauSamples);
% 
%     % 初始化模糊函数矩阵 AF: 大小 Ntau x length(vAxis)
%     AF = zeros(Ntau, length(vAxis));
% 
%     % 遍历速度
%     for iv = 1:length(vAxis)
%         v = vAxis(iv);
%         alpha = 1 + (v/c); % 多普勒伸缩因子
% 
%         % 计算 p, q 用于对信号进行重采样
%         [p, q] = rat(alpha, 1e-7); 
%         % 对 x 进行重采样近似多普勒效应
%         xDoppler = resample(x, p, q); 
%         ND = length(xDoppler);
% 
%         % 对每个延迟 tau 进行相关运算
%         for it = 1:Ntau
%             tau = tauSamples(it);
%             if tau >= 0
%                 lenC = min(ND - tau, Nx);
%                 if lenC > 0
%                     segDopp = xDoppler(1+tau : tau+lenC);
%                     segX    = x(1 : lenC);
%                     AF(it, iv) = abs(segDopp * segX');
%                 else
%                     AF(it, iv) = 0;
%                 end
%             else
%                 tauPos = -tau;
%                 lenC = min(ND, Nx - tauPos);
%                 if lenC > 0
%                     segDopp = xDoppler(1 : lenC);
%                     segX    = x(1+tauPos : tauPos+lenC);
%                     AF(it, iv) = abs(segDopp * segX');
%                 else
%                     AF(it, iv) = 0;
%                 end
%             end
%         end
%     end
% 
%     % 将 AF 最大值归一化为 1（可选）
%     AF = AF ./ max(AF(:));
% 
%     % 将延迟索引转换为时间 (秒)
%     tauAxis = tauSamples / fs;
% end
% 宽带模糊函数计算函数
function AF = widebandAF(s1, s2, tau_range, d_range, Fs)
    N1 = length(s1);
    N2 = length(s2);
    t1 = (0:N1-1)/Fs;
    t2 = (0:N2-1)/Fs;
    
    AF = zeros(length(tau_range), length(d_range));
    
    parfor i_tau = 1:length(tau_range)
        tau = tau_range(i_tau);
        AF_row = zeros(1, length(d_range));
        for i_d = 1:length(d_range)
            d = d_range(i_d);
            interp_points = d*(t1 - tau);
            s2_transformed = interp1(t2, s2, interp_points, 'linear', 0);
            corr = xcorr(s1(:), s2_transformed(:));
            AF_row(i_d) = max(abs(corr));
        end
        AF(i_tau, :) = AF_row;
    end
end