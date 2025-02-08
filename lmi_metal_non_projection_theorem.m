%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%% 滤波器参数求解 模态无关李雅普诺夫函数 不使用投影定理  %%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
%% 实值矩阵维数参数（可调节）
nw = 1; % 扰动维数
nz = 1; % 待估计信号维数
ny = 1; % 测量信号维数 ny = nz
nhx = 2;  nvx = 3;   % 系统维数（水平和垂直）
nhf = 2;  nvf = 3;   % 滤波器维数（水平和垂直）
nx = nhx + nvx; % 系统总维数
nf = nhf + nvf; % 滤波器总维数
%% 数据来源：H∞ Model Reduction for 2-D Discrete Markovian Jump Systems --Khalid Badie
% 模态实值矩阵 模态1
A1 = [0.3846   -0.1479    0.2959   -0.1479    0.0740;
           0         0    1.0000         0         0;
      0.3846   -0.1479    0.2959   -0.1479    0.0740;
           0         0    1.0000         0         0;
           0    1.0000         0         0         0];    % nx nx
B1 = [0.9 1.2 0.8 0.17 0.4]';  % nx nw
C1 = [0 0 1 0 0];    % ny nx
D1 = 1;   % ny nw
J1 = [0 1 1 1 0]; % nz nx
L1 = 0.8;   % nz nw
% 模态实值矩阵 模态2
A2 = [0.4320   -0.1479    0.2959   -0.1479    0.0740;
           0         0    1.0000         0         0;
      0.4320   -0.1479    0.2959   -0.1479    0.0740;
           0         0    1.0000         0         0;
           0    1.0000         0         0         0];
B2 = [-0.6 0.7 1.3 1.1 -1.4]';
C2 = C1;  D2 = D1;  J2 = J1;  L2 = L1;
%% 系统伯努利过程 alpha 系统切换
Phi = [0.4   0.6];  a1 = 0.4;  a2 = 0.6;
%% 条件概率过程 beta 
beta = [0.8  0.2; 0.3  0.7]; % （异步）异步滤波器切换
% beta = eye(2,2); % （同步）滤波器存在切换 追随系统模态
% beta = [1 0; 1 0]; % （模式无关）滤波器不存在切换 仅存在模态一
%% Sensor failures 传感器故障模型（随着阶数变化）
% 传感器部分故障  0.2 < eta < 0.4  
% etaUp = 0.3;  etaDown = 0.1;
% 传感器部分故障  0.4 < eta < 0.6  
% etaUp = 0.5;  etaDown = 0.1;
% 传感器部分故障  0.6 < eta < 0.8  
etaUp = 0.7;  etaDown = 0.1;
% 传感器部分故障  0.8 < eta < 1  
% etaUp = 0.9;  etaDown = 0.1;
% 传感器无故障    
% etaUp = 1;  etaDown = 0;
%% 决策变量
gamma = sdpvar(1, 1, 'symmetric');  %性能指标
epsilon = sdpvar(1, 1, 'full');  % 放缩参数1
tau = sdpvar(1, 1, 'full');  % 放缩参数2
trans = [eye(nf) zeros(nf, nx-nf)]'; % 降阶滤波调节矩阵Tao8

Ph1 = sdpvar(nhx, nhx, 'symmetric'); % sum(1:nhx)
Pv1 = sdpvar(nvx, nvx, 'symmetric'); % sum(1:nvx)
Ph2 = sdpvar(nhf, nhf, 'symmetric'); % sum(1:nhf)
Pv2 = sdpvar(nvf, nvf, 'symmetric'); % sum(1:nvf)
Ph3 = sdpvar(nhx, nhf, 'full'); % nhx*nhf
Pv3 = sdpvar(nvx, nvf, 'full'); % nvx*nvf

P1 = blkdiag(Ph1, Pv1);  P2 = blkdiag(Ph2, Pv2);  P3 = blkdiag(Ph3, Pv3);
P = [P1  P3;
     P3' P2];
%% 系统模态/滤波器模态 循环定义决策变量
NumSysMode = 2;  NumFilterMode = 2;  % 系统模态个数 滤波器模态个数
for m = 1 : NumSysMode
    for n = 1 : NumFilterMode
       %% 额外决策变量
        eval(['Vmn_1_', int2str(m), int2str(n), '= sdpvar(nx, nx, ''symmetric'');']); % sum(1:nx)
        eval(['Vmn_2_', int2str(m), int2str(n), '= sdpvar(nf, nx, ''full'');']); % nf*nx
       %% 决策变量
        eval(['Rmn_', int2str(m), int2str(n), '= sdpvar(nw, nw, ''symmetric'');']); % sum(1:nw)
        eval(['Gmn_', int2str(m), int2str(n), '= sdpvar(nx+nf, nx+nf, ''symmetric'');']); % sum(1:nx+nf)
        eval(['Afn_', int2str(n), '= sdpvar(nf, nf, ''full'');']);  % nf^2
        eval(['Bfn_', int2str(n), '= sdpvar(nf, ny, ''full'');']);  % nf*ny
        eval(['Jfn_', int2str(n), '= sdpvar(nz, nf, ''full'');']);  % nz*nf
        eval(['Lfn_', int2str(n), '= sdpvar(nz, ny, ''full'');']);  % nz*ny
        eval(['Vn_', int2str(n), '= sdpvar(nf, nf, ''full'');']);  % nf^2 
    end
end
%% 循环定义LMI
for m = 1:NumSysMode
    for n = 1:NumFilterMode
       %% 决策变量
        Am = eval(['A', int2str(m)]);  Bm = eval(['B', int2str(m)]);  Cm = eval(['C', int2str(m)]);
        Dm = eval(['D', int2str(m)]);  Jm = eval(['J', int2str(m)]);  Lm = eval(['L', int2str(m)]); % 系统参数
        Afn = eval(['Afn_', int2str(n)]);  Bfn = eval(['Bfn_', int2str(n)]);  
        Jfn = eval(['Jfn_', int2str(n)]);  Lfn = eval(['Lfn_', int2str(n)]); % 参数矩阵
        Vn = eval(['Vn_', int2str(n)]); % 解耦矩阵
        Rmn = eval(['Rmn_', int2str(m), int2str(n)]);  Gmn = eval(['Gmn_', int2str(m), int2str(n)]); % 松弛矩阵
        Vmn_1 = eval(['Vmn_1_', int2str(m), int2str(n)]);  Vmn_2 = eval(['Vmn_2_', int2str(m), int2str(n)]);  
        bar_Vmn = [Vmn_1;  Vmn_2]; % bar_Vmn
        Vmn = [Vmn_1  trans*Vn;
               Vmn_2        Vn]; % Vmn
       %% 论文中的（21）
        Jmn = [Jm-Lfn*etaUp*Cm  -Jfn]; % Jmn 
        hat_Cm = [Cm zeros(ny, nf)]; % hat_Cm 
        Sigma_1 = [-eye(nz)    Jmn;
                    Jmn'      -Gmn]; % Sigma_1 
        Sigma_2 = [zeros(nx,nz)  Vmn_1*Am + trans*Bfn*etaUp*Cm   trans*Afn;
                   zeros(nf,nz)  Vmn_2*Am +       Bfn*etaUp*Cm         Afn]; % Sigma_2            
        Sigma_3 = [zeros(nx,ny)   trans*Bfn;
                   zeros(nf,ny)         Bfn]; % Sigma_3      
        Sigma_4 = [zeros(nz,ny)        -Lfn;
                   tau*etaDown*hat_Cm'  zeros(nx+nf,ny)]; % Sigma_4      
        Sigma = [P - Vmn - Vmn'  Sigma_2    Sigma_3;
                 Sigma_2'        Sigma_1    Sigma_4;
                 Sigma_3'        Sigma_4'  -tau*eye(2*ny)]; % Sigma  
       %% 论文中的（22）
        Upsilon_1 = [-eye(nz)             (Lm-Lfn*etaUp*Dm);
                     (Lm-Lfn*etaUp*Dm)'  -Rmn]; % Omega1_st 
        Upsilon_2 = [zeros(nx,nz)  Vmn_1*Bm + trans*Bfn*etaUp*Dm;
                     zeros(nf,nz)  Vmn_2*Bm +       Bfn*etaUp*Dm];  % Omega2_st         
        Upsilon_3 = [zeros(nx,ny)  trans*Bfn;
                     zeros(nf,ny)        Bfn]; % Omega3_st          
        Upsilon_4 = [zeros(nz,ny)        -Lfn;
                     epsilon*etaDown*Dm'  zeros(nw,ny)]; % Omega4_st    
        Upsilon = [P- Vmn - Vmn'  Upsilon_2    Upsilon_3;
                   Upsilon_2'     Upsilon_1    Upsilon_4;
                   Upsilon_3'     Upsilon_4'  -epsilon*eye(2*ny)]; % Omega st     
        eval(['Condition_', int2str(m), int2str(n),'= [Sigma <= 0, Upsilon <= 0, Gmn >= 0, Rmn >= 0, P >= 0];']);
    end
end
%% LMIs 论文中的（10） （11）
fia11 = a1*beta(1,1);  fia12 = a1*beta(1,2);  fia21 = a2*beta(2,1);  fia22 = a2*beta(2,2); 
TempRmn = zeros(nw); TempGmn = zeros(nx+nf);
for m = 1 : NumSysMode
    for n = 1 : NumFilterMode
        fia = eval(['fia', int2str(m), int2str(n)]);
        Rmn = eval(['Rmn_', int2str(m), int2str(n)]);
        Gmn = eval(['Gmn_', int2str(m), int2str(n)]);
        TempRmn = TempRmn + fia*trace(Rmn);
        TempGmn = TempGmn + fia*Gmn;
    end
end
LMIRmn = TempRmn - gamma;
LMIGmn = TempGmn - P;
%% 线性约束条件，需要同时满足
Constraint = [LMIRmn <= 0, LMIGmn <= 0, gamma >= 0, epsilon >= 0, tau >= 0;
              Condition_11, Condition_12, Condition_21, Condition_22];
sdpsettings('solver', 'mosek', 'verbos', 0); % 设置求解器为mosek，并打印少量信息
reuslt = optimize(Constraint, gamma); % 自动求解的优化函数
if reuslt.problem == 0 %problem = 0 代表求解成功
    [primal,~]=check(Constraint); % 检查约束
    if min(primal)>=0 && all(primal([1:2,4])>0) % 判断残差值
        disp('***********************************');
        disp('****Constraints are guaranteed*****');
        disp('**Primal residual is/are positive**');
        disp('***********************************');
        Gamma = sqrt(value(gamma))
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        Gamma = sqrt(value(gamma))
        check(Constraint); % 检查残余值
    end
else
    disp('**************************************');
    disp('****Constraints are not guaranteed****');
    disp('**************************************'); 
    reuslt.info
    yalmiperror(reuslt.problem) % 打印出错信息
end
disp(['运行时间: ', num2str(toc)]); 