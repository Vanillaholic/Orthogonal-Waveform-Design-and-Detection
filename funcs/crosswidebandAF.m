function [CAF, tauAxis, vAxis] = crosswidebandAF(x, y, fs, c,vAxis,maxDelay)
%CROSSWIDEBANDAF 计算两个信号的宽带互模糊函数
%
%   [CAF, tauAxis, vAxis] = crosswidebandAF(x, y, fs, vAxis, c, maxDelay)
%
% 输入参数：
%   x       : 参考信号（行向量或列向量）
%   y       : 待检测信号（行向量或列向量）
%   fs      : 采样率（Hz）
%   vAxis   : 速度向量，例如：-20:1:20（单位 m/s）
%   c       : 波速（例如声速 1500 m/s 或光速 3e8 m/s）
%   maxDelay: 最大延迟（以采样点计），默认为 length(x)-1
%
% 输出参数：
%   CAF     : 互模糊函数矩阵，行对应延迟（tau），列对应速度
%   tauAxis : 延迟轴（单位：秒），由采样点转换而来
%   vAxis   : 与输入相同的速度向量
%
% 原理说明：
%   对于每个速度 v，计算多普勒伸缩因子 alpha = 1 + v/c，
%   利用 rat() 将 alpha 用分数 (p,q) 表示，然后使用 resample()
%   对信号 y 进行重采样以模拟多普勒效应，得到 y_Doppler。
%   接下来计算 x 与 y_Doppler 在不同延迟下的互相关，
%   取其幅值作为互模糊函数的值。
%

    % 确保信号为行向量
    x = x(:).'; 
    y = y(:).';

    if nargin < 6
        maxDelay = length(x) - 1;
    end

    % 延迟采样索引，从 -maxDelay 到 maxDelay
    tauSamples = -maxDelay:maxDelay;
    Ntau = length(tauSamples);
    Nv = length(vAxis);
    
    CAF = zeros(Ntau, Nv); % 初始化互模糊函数矩阵

    % 对每个速度计算互模糊函数
    for iv = 1:Nv
        v = vAxis(iv);
        alpha = 1 + (v / c);   % 多普勒伸缩因子
        % 用 rat 将 alpha 表示为 p/q 分数（可调逼近精度）
        [p, q] = rat(alpha, 1e-7);
        % 重采样信号 y，模拟多普勒效应：当 alpha >1时扩展，当 alpha <1时压缩
        yDoppler = resample(y, p, q);
        
        % 计算 x 与 yDoppler 的互相关，范围为 -maxDelay ~ maxDelay
        ccf = xcorr(x, yDoppler, maxDelay);
        % 存储互相关幅值
        CAF(:, iv) = abs(ccf(:));
    end

    % 归一化（可选），将最大值归一化为1
    CAF = CAF / max(CAF(:));
    
    % 将延迟采样索引转换为秒
    tauAxis = tauSamples / fs;
end
