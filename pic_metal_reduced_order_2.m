%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% 降阶滤波器状态轨迹图 nhf=1 nhv=1 %%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
%% 绘图坐标 绘图线条 坐标轴字体 大小
fontSizeXY = 28;  lineWidth = 1;  fontSizeAxis = 23; markerSize = 20;
T = 40;  dt = 1;  Nt = T/dt; % 参数：水平方向
Lx = 40;  dx = 1;  Nx = Lx/dx; % 参数：垂直方向
rand('state', 50); % 随机数生成器 rand 状态保持不变
%% 系统参数来源：H∞ Model Reduction for 2-D Discrete Markovian Jump Systems --Khalid Badie
% 模态实值矩阵 模态1
As1 = [0.3846   -0.1479    0.2959   -0.1479    0.0740;
            0         0    1.0000         0         0;
       0.3846   -0.1479    0.2959   -0.1479    0.0740;
            0         0    1.0000         0         0;
            0    1.0000         0         0         0];    % nx nx
Bs1 = [0.9 1.2 0.8 0.17 0.4]';  % nx nw
Cs1 = [0 0 1 0 0];    % ny nx
Ds1 = 1;   % ny nw
Js1 = [0 1 1 1 0]; % nz nx
Ls1 = 0.8;   % nz nw
% 模态实值矩阵 模态2
As2 = [0.4320   -0.1479    0.2959   -0.1479    0.0740;
            0         0    1.0000         0         0;
       0.4320   -0.1479    0.2959   -0.1479    0.0740;
            0         0    1.0000         0         0;
            0    1.0000         0         0         0];
Bs2 = [-0.6 0.7 1.3 1.1 -1.4]';
Cs2 = Cs1;  Ds2 = Ds1;  Js2 = Js1;  Ls2 = Ls1;
%% (异步) H2 异步双模滤波器参数
% 异步滤波器模态1
At1 = [-0.0338    0.0270;
       -0.0097    0.0077];
Bt1 = [-3.9793;
       -1.1356];
Jt1 = [-0.2325   -0.8841];
Lt1 = 1.2136;
% 异步滤波器模态2
At2 = [ 0.0713    0.0228;
        0.0203    0.0065];
Bt2 = [ -3.9620;
        -1.1305];
Jt2 = [-0.2435   -0.8830];
Lt2 = 1.2366;
%% 系统伯努利过程 alpha 系统切换
Phi = [0.4   0.6];  a1 = 0.4;  a2 = 0.6;
%% （异步） 条件概率过程 beta 异步滤波器切换
beta = [0.8  0.2; 0.3  0.7];
%% 仿真迭代存储矩阵
eta = 0.7; % 传感器故障强度
x1 = zeros(Nt, Nx);  x2 = zeros(Nt, Nx);  x3 = zeros(Nt, Nx);  x4 = zeros(Nt, Nx);  x5 = zeros(Nt, Nx); % 系统输出序列存储数组(滤波前)
xf1 = zeros(Nt, Nx);  xf2 = zeros(Nt, Nx); % 滤波输出序列存储数组(滤波后)
y = zeros(Nt, Nx); % 系统测量输出序列存储数组
yf = zeros(Nt, Nx); % 系统故障测量输出序列存储数组
z = zeros(Nt, Nx); % 待估计的信号(实际值)序列存储数组
zsim = zeros(Nt, Nx); % 待估计的信号(估计值)序列存储数组
zerror = zeros(Nt, Nx); % 估计误差序列存储数组
zsum = zeros(Nt, Nx); % H2 范数序列存储数组
w = zeros(Nt,Nx); w(1,1) = 1; % 冲激信号 impulse signal
Sysseq = zeros(Nt+1,Nx+1); % 系统状态序列存储数组
Filterseq = zeros(Nt+1,Nx+1); % 异步滤波器状态序列存储数组
%% 系统 滤波器 边界条件
% x1(1, 1:Nx) = 3.3;  x2(1, 1:Nx) = 2.7;  xf1(1, 1:Nx) = -2.5; % 水平方向
% x3(1:Nt, 1) = -1.5;  x4(1:Nt, 1) = -2.0;  x5(1:Nt, 1) = -2.1;  xf2(1:Nt, 1) = 1.8; % 垂直方向
%% 系统 滤波器 模态初值 边界条件
Sysseq(:,1) = 1;  Sysseq(1,:) = 2; % 系统模态初值
Filterseq(:,1) = 1;  Filterseq(1,:) = 2; % 滤波器模态初值
%% 系统与异步滤波器模态切换代码
for i = 1:Nt+1
    for j = 1:Nx+1
            a = rand;  b = rand;
                if a < a1        
                    flagSys = 1;       % System jump to mode 1
                    if b < beta(1,1)
                        flagFilter = 1;
                        else
                        flagFilter = 2;
                    end
                else      
                    flagSys = 2;       % System jump to mode 2
                        if b < beta(2,1)
                        flagFilter = 1;
                        else
                        flagFilter = 2;
                        end
                end
        Sysseq(i,j) = flagSys;  Filterseq(i,j) = flagFilter;
    end
end
%% 2D MJSs 动态方程
for i = 1:Nt
    for j = 1:Nx
        switch(Sysseq(i,j)) % 系统模式切换
            case 1
                As = As1;  Bs = Bs1;  Cs = Cs1;  Ds = Ds1;  Js = Js1;  Ls = Ls1;
            case 2
                As = As2;  Bs = Bs2;  Cs = Cs2;  Ds = Ds2;  Js = Js2;  Ls = Ls2;
        end
        switch(Filterseq(i,j)) % 双模异步滤波器参数切换
            case 1
                At = At1;  Bt = Bt1;  Jt = Jt1;  Lt = Lt1;
        	case 2
                At = At2;  Bt = Bt2;  Jt = Jt2;  Lt = Lt2;
        end       
        x1(i+1,j) = As(1,:)*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Bs(1,1)*w(i,j);
        x2(i+1,j) = As(2,:)*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Bs(2,1)*w(i,j);
        x3(i,j+1) = As(3,:)*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Bs(3,1)*w(i,j);
        x4(i,j+1) = As(4,:)*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Bs(4,1)*w(i,j);
        x5(i,j+1) = As(5,:)*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Bs(5,1)*w(i,j);  % 系统状态 五向量
        y(i,j) = Cs*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Ds*w(i,j); % 输出信号  
        z(i,j) = Js*[x1(i,j); x2(i,j); x3(i,j); x4(i,j); x5(i,j)] + Ls*w(i,j); % 待测信号    
        yf(i,j) =  eta*y(i,j); % 故障输出信号    
        xf1(i+1,j) = At(1,:)*[xf1(i,j); xf2(i,j)] + Bt(1,1)*yf(i,j);
        xf2(i,j+1) = At(2,:)*[xf1(i,j); xf2(i,j)] + Bt(2,1)*yf(i,j); % 全阶滤波器信号 
        zsim(i,j) = Jt*[xf1(i,j); xf2(i,j)] + Lt*yf(i,j); % 估计信号
        zerror(i,j) = z(i,j) - zsim(i,j); % 滤波误差
        zsum(i,j) = sqrt(norm(zerror(1:i, 1:j), 2)); % 误差信号的2范数
    end
end
%% 滤波误差
figure('name', '2 阶 滤波误差信号');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = zerror(1:dt:T, 1:dt:Lx);  surf(x, y, z);  
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$i$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
zlabel('$\tilde{z}_{i,j}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
axis([0 T 0 Lx -2 4]);  grid on;  grid minor;  box on;
%% 性能比较
figure('name', '2 阶 性能比较');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = zsum(1:dt:T, 1:dt:Lx);  surf(x, y, z);
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$q$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$p$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
zlabel('$\gamma_{p,q}$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
axis([0 T 0 Lx 0.2 1]);  grid on;  grid minor;  box on;
%% H2 范数求解
disp('****************************');
disp('****水平/垂直 爱趣 2 范数****');
disp('****************************');
Norm_order_2 = sqrt(norm(zerror(1:i,1:j), 2))
disp(['运行时间: ', num2str(toc)]); 