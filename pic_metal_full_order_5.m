%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% 全阶滤波器状态轨迹图 nhf=2 nhv=3 %%%%%%%%%%%%%%%%%%%%%
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
At1 = [ 0.1082   -0.0957    0.0767   -0.2368    0.1268;
        0.0006   -0.0020    0.0055    0.0008   -0.0004;
        0.0765   -0.0580   -0.4141   -0.2681    0.1346;
        0.0086   -0.0227    0.3534    0.0091   -0.0047;
       -0.0155    0.9856    0.2030    0.0018    0.0007];
Bt1 = [-0.4180;
       -1.4160;
       -1.2235;
       -0.8718;
        0.3259];
Jt1 = [0.0214   -1.0086   -0.1511   -0.9911   -0.0042];
Lt1 = 1.2328;
% 异步滤波器模态2
At2 = [ 0.1074   -0.0549    0.6001   -0.2560    0.1375;
        0.0047   -0.0074    0.1492    0.0035   -0.0018;
        0.0730   -0.0519   -0.5524   -0.2783    0.1385;
        0.0016   -0.0061    0.0910    0.0028   -0.0008;
       -0.0020    0.9816    0.7468    0.0068   -0.0013];
Bt2 = [ 0.2325;
       -1.1985;
       -1.4352;
       -1.2844;
        1.1114];
Jt2 = [0.0220   -1.0090   -0.1519   -0.9909   -0.0044];
Lt2 = 1.2324;
%% 系统伯努利过程 alpha 系统切换
Phi = [0.4   0.6];  a1 = 0.4;  a2 = 0.6;
%% （异步） 条件概率过程 beta 异步滤波器切换
beta = [0.8  0.2; 0.3  0.7];
%% 仿真迭代存储矩阵
eta = 0.7; % 传感器故障强度
x1 = zeros(Nt, Nx);  x2 = zeros(Nt, Nx);  x3 = zeros(Nt, Nx);  x4 = zeros(Nt, Nx);  x5 = zeros(Nt, Nx); % 系统输出序列存储数组(滤波前)
xf1 = zeros(Nt, Nx);  xf2 = zeros(Nt, Nx);  xf3 = zeros(Nt, Nx);  xf4 = zeros(Nt, Nx);  xf5 = zeros(Nt, Nx); % 滤波输出序列存储数组(滤波后)
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
x1(1, 1:Nx) = 3.3;  x2(1, 1:Nx) = 2.7;  xf1(1, 1:Nx) = -2.5;  xf2(1, 1:Nx) = -2.2; % 水平方向
x3(1:Nt, 1) = -1.5;  x4(1:Nt, 1) = -2.0;  x5(1:Nt, 1) = -2.1;  xf3(1:Nt, 1) = 1.8;  xf4(1:Nt, 1) = 1.7;  xf5(1:Nt, 1) = 2.3; % 垂直方向
%% 系统 滤波器 模态初值 边界条件
Sysseq(:,1) = 1;  Sysseq(1,:) = 2; % 系统模态初值
Filterseq(:,1) = 1;  Filterseq(1,:) = 2; % 滤波器模态初值
%% 系统与异步滤波器模态切换代码
for i = 2:Nt+1
    for j = 2:Nx+1
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
%% 2D BJSs 动态方程迭代过程
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
        yf(i,j) = eta*y(i,j); % 故障输出信号
        xf1(i+1,j) = At(1,:)*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); x5(i,j)] + Bt(1,1)*yf(i,j);
        xf2(i+1,j) = At(2,:)*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); x5(i,j)] + Bt(2,1)*yf(i,j);
        xf3(i,j+1) = At(3,:)*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); x5(i,j)] + Bt(3,1)*yf(i,j);
        xf4(i,j+1) = At(4,:)*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); x5(i,j)] + Bt(4,1)*yf(i,j);
        xf5(i,j+1) = At(5,:)*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); x5(i,j)] + Bt(5,1)*yf(i,j); % 全阶滤波器信号
        zsim(i,j) = Jt*[xf1(i,j); xf2(i,j); xf3(i,j); xf4(i,j); xf5(i,j)] + Lt*yf(i,j); % 估计信号
        zerror(i,j) = z(i,j) - zsim(i,j); % 滤波误差
        zsum(i,j) = sqrt(norm(zerror(1:i, 1:j), 2)); % 误差信号的2范数
    end
end
%% 待测信号
% figure('name', '实际信号');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = z(1:dt:T,1:dt:Lx);  surf(x,y,z);
% xlabel('$j$','Interpreter','latex','Fontsize',16);
% ylabel('$i$','Interpreter','latex','Fontsize',16);
% zlabel('${z}_{i,j}$','Interpreter','latex','Fontsize',16');
% % axis([0 T 0 Nx -5 10]);  grid on;  grid minor;  box on;
%% 估计信号
% figure('name', '估计信号');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = zsim(1:dt:T,1:dt:Lx);  surf(x,y,z);
% xlabel('$j$','Interpreter','latex','Fontsize',16);
% ylabel('$i$','Interpreter','latex','Fontsize',16);
% zlabel('$\hat{z}_{i,j}$','Interpreter','latex','Fontsize',16');
% % axis([0 T 0 Nx -5 10]);  grid on;  grid minor;  box on;
%% 滤波误差
figure('name', '5 阶 滤波误差信号');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = zerror(1:dt:T, 1:dt:Lx);
surf(x, y, z);
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$i$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
zlabel('$\tilde{z}_{i,j}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
axis([0 T 0 Lx -2 4]);  grid on;  grid minor;  box on;
%% 性能比较
figure('name', '5 阶 性能比较');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = zsum(1:dt:T, 1:dt:Lx);  surf(x, y, z);
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$q$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$p$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
zlabel('$\gamma_{p,q}$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
axis([0 T 0 Lx 0.2 0.82]);  grid on;  grid minor;  box on;
%% 系统模式切换
figure('name', '系统模式切换信号');
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);  z = Sysseq(1:dt:T+1, 1:dt:Lx+1);
surf(x, y, z);  view(0, 90);  axis([0 T 0 Lx 0.5 4.5]);
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$i$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
%% 滤波器模式切换
figure('name', '滤波器模式切换信号');
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);  z = Filterseq(1:dt:T+1, 1:dt:Lx+1);
surf(x, y, z);  view(0, 90);  axis([0 T 0 Lx 0.5 4.5]);
set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
xlabel('$j$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
ylabel('$i$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
%% H2 范数求解
disp('****************************');
disp('****水平/垂直 爱趣 2 范数****');
disp('****************************');
Norm_order_5 = sqrt(norm(zerror(1:i,1:j), 2))
disp(['运行时间: ', num2str(toc)]); 