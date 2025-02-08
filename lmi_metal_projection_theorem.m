%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% 滤波器参数求解 模态无关李雅普诺夫函数 使用投影定理 %%%%%%%%%%%%%%%
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
% gamma = sdpvar(1, 1, 'symmetric');  % 性能指标
gamma = 1.0377^2;
epsilon = sdpvar(1, 1, 'symmetric');  % 放缩参数1
tau = sdpvar(1, 1, 'symmetric');  % 放缩参数2
Std = [eye(nf)  zeros(nf, nx-nf)]'; % 降阶滤波调节矩阵Tao

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
for m = 1:NumSysMode
    for n = 1:NumFilterMode
       %% 额外决策变量
        eval(['Vmn_1_', int2str(m), int2str(n), '= sdpvar(nx, nx, ''symmetric'');']); % 验证正交补充相乘为 0 不计入 LMIs
        eval(['Vmn_2_', int2str(m), int2str(n), '= sdpvar(nf, nx, ''full'');']); % 验证正交补充相乘为 0 不计入 LMIs
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
for m = 1 : NumSysMode
    for n = 1 : NumFilterMode
        Vmn_1 = eval(['Vmn_1_', int2str(m), int2str(n)]);  Vmn_2 = eval(['Vmn_2_', int2str(m), int2str(n)]);  
        bar_Vmn = [Vmn_1;  Vmn_2]; % bar_Vmn
       %% 决策变量
        Am = eval(['A', int2str(m)]);  Bm = eval(['B', int2str(m)]);  Cm = eval(['C', int2str(m)]);
        Dm = eval(['D', int2str(m)]);  Jm = eval(['J', int2str(m)]);  Lm = eval(['L', int2str(m)]); % 系统参数
        Afn = eval(['Afn_', int2str(n)]);  Bfn = eval(['Bfn_', int2str(n)]);  
        Jfn = eval(['Jfn_', int2str(n)]);  Lfn = eval(['Lfn_', int2str(n)]); % 参数矩阵
        Vn = eval(['Vn_', int2str(n)]);
        Rmn = eval(['Rmn_', int2str(m), int2str(n)]);  Gmn = eval(['Gmn_', int2str(m), int2str(n)]); % 松弛矩阵
        Vmn = [zeros(nx, nx)  Std*Vn;
               zeros(nf, nx)      Vn]; % Vmn
       %% 论文中的 Sigma_delta
        Jmn = [Jm - Lfn*etaUp*Cm  -Jfn];  % Jmn
        hat_Cm = [Cm zeros(ny, nf)]; % hat_Cm 
        Sigma_1 = [-eye(nz)    Jmn;
                    Jmn'      -Gmn]; % Sigma_1    
        Sigma_2 = [zeros(nx,nz)  Std*Bfn*etaUp*Cm  Std*Afn;
                   zeros(nf,nz)      Bfn*etaUp*Cm      Afn]; % Sigma_2                
        Sigma_3 = [zeros(nx,ny)  Std*Bfn;
                   zeros(nf,ny)      Bfn]; % Sigma_3
        Sigma_4 = [zeros(nz,ny)         -Lfn;
                   tau*etaDown*hat_Cm'   zeros(nx+nf,ny)]; % Sigma_4    
        Sigma_delta = [P - Vmn - Vmn'  Sigma_2    Sigma_3;
                       Sigma_2'        Sigma_1    Sigma_4;
                       Sigma_3'        Sigma_4'  -tau*eye(2*ny)]; % Sigma_delta 
       %% 论文中的 Upsilon_delta
        Upsilon_1 = [-eye(nz)            (Lm-Lfn*etaUp*Dm);
                     (Lm-Lfn*etaUp*Dm)'  -Rmn]; % Upsilon_1 
        Upsilon_2 = [zeros(nx,nz)  Std*Bfn*etaUp*Dm;
                     zeros(nf,nz)      Bfn*etaUp*Dm]; % Upsilon_2        
        Upsilon_3 = [zeros(nx,ny)  Std*Bfn;
                     zeros(nf,ny)      Bfn]; % Upsilon_3     
        Upsilon_4 = [zeros(nz,ny)          -Lfn;
                     epsilon*etaDown*Dm'    zeros(nw,ny)]; % Upsilon_4     
        Upsilon_delta = [P - Vmn - Vmn'  Upsilon_2    Upsilon_3;
                         Upsilon_2'      Upsilon_1    Upsilon_4;
                         Upsilon_3'      Upsilon_4'  -epsilon*eye(2*ny)]; % Upsilon_delta     
       %% ker_mathscr_A 和 ker_mathscr_B 分别为 mathscr_A 和 mathscr_B 的正交补
        ASims = [zeros(nx, nf)  zeros(nx, nz)  Am  zeros(nx, nf)  zeros(nx, 2*ny)]; % 2*ny 是因为引入故障传感器而产生的
        ker_mathscr_A = [ASims;  eye(nx+2*nf+nz+2*ny)];  % 可以验证 mathscr_A*ker_mathscr_A = 0
        BSims = [zeros(nx, nf)  zeros(nx, nz)  Bm  zeros(nx, 2*ny)]; % 2*ny 是因为引入故障传感器而产生的
        ker_mathscr_B = [BSims;  eye(nf+nz+2*ny+nw)];    % 可以验证 mathscr_B*ker_mathscr_B = 0
       %% ker_mathscr_V 和 ker_mathscr_U 分别是 mathscr_V 和 mathscr_U 的正交补
        ker_mathscr_V = [zeros(nx+nf, nx+nf+nz+2*ny);  eye(nx+nf+nz+2*ny)];   % 可以验证 mathscr_V*ker_mathscr_V = 0
        ker_mathscr_U = [zeros(nx+nf, nz+nw+2*ny);  eye(nz+nw+2*ny)];    % 可以验证 mathscr_U*ker_mathscr_U = 0
       %% 论文中的 (26) (27) (28) (29)
        LMI_26 = ker_mathscr_A'*Sigma_delta*ker_mathscr_A; % Theorem 3 (26)
        LMI_27 = ker_mathscr_V'*Sigma_delta*ker_mathscr_V; % Theorem 3 (27)
        LMI_28 = ker_mathscr_B'*Upsilon_delta*ker_mathscr_B; % Theorem 3 (28)
        LMI_29 = ker_mathscr_U'*Upsilon_delta*ker_mathscr_U; % Theorem 3 (29)
       %% ★★★验证★★★：用于验证正交补充相乘为 0 未使用投影定理
        mathscr_A = [-eye(nx)  zeros(nx,nf)  zeros(nx,nz)  Am  zeros(nx,nf)  zeros(nx,2*ny)]; % 2*ny 是因为引入故障传感器而产生的
        mathscr_V = [eye(nx+nf) zeros(nx+nf, nx+nf+nz+2*ny)];  
        mathscr_B = [-eye(nx)  zeros(nx,nf)  zeros(nx,nz)  Bm  zeros(nx,2*ny)]; % 2*ny 是因为引入故障传感器而产生的
        mathscr_U = [eye(nx+nf) zeros(nx+nf, nz+nw+2*ny)];     
        LMI_32 = Sigma_delta + mathscr_A'*bar_Vmn'*mathscr_V + mathscr_V'*bar_Vmn*mathscr_A; % (32)
        LMI_33 = Upsilon_delta + mathscr_B'*bar_Vmn'*mathscr_U + mathscr_U'*bar_Vmn*mathscr_B; % (33)
        eval(['Condition_', int2str(m), int2str(n), '= [LMI_26 <= 0, LMI_27 <= 0, LMI_28 <= 0, LMI_29 <= 0, Gmn >= 0, Rmn >= 0, P>=0];']);
    end
end
%% LMIs 论文中的（10） （11）
fia11 = a1*beta(1,1);  fia12 = a1*beta(1,2);  fia21 = a2*beta(2,1);  fia22 = a2*beta(2,2);
TempRmn = zeros(nw);  TempGmn = zeros(nx+nf);
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
        Decision_Variables = sum(1:nhx) + sum(1:nvx) + sum(1:nhf) + sum(1:nvf) + nhx*nhf + nvx*nvf...
                             + NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx + sum(1:nw) + sum(1:nx+nf))...
                             + NumFilterMode*(nf^2 + nf*ny + nz*nf + nz*ny + nf^2) + 2
        Extra_Variables = NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx)
        Decision_Variables_with_Projection_Theorem = Decision_Variables - NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx)
        Reduction_Rate = sprintf('%.1f%%',(NumSysMode*NumSysMode*(sum(1:nx)+nf*nx)/Decision_Variables*100))
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        Gamma = sqrt(value(gamma))
        Decision_Variables = sum(1:nhx) + sum(1:nvx) + sum(1:nhf) + sum(1:nvf) + nhx*nhf + nvx*nvf...
                             + NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx + sum(1:nw) + sum(1:nx+nf))...
                             + NumFilterMode*(nf^2 + nf*ny + nz*nf + nz*ny + nf^2) + 2
        Extra_Variables = NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx)
        Decision_Variables_with_Projection_Theorem = Decision_Variables - NumSysMode*NumFilterMode*(sum(1:nx) + nf*nx)
        Reduction_Rate = sprintf('%.1f%%',(NumSysMode*NumSysMode*(sum(1:nx)+nf*nx)/Decision_Variables*100))
%         check(Constraint); % 检查残余值
    end
else
    disp('**************************************');
    disp('****Constraints are not guaranteed****');
    disp('**************************************'); 
    reuslt.info
    yalmiperror(reuslt.problem) % 打印出错信息
end
%% 滤波器参数矩阵
disp('****************************');
disp('****滤波器参数矩阵 模态 1****');
disp('****************************');
Af1 = value(inv(value(Vn_1))*Afn_1), Bf1 = value(inv(value(Vn_1))*Bfn_1), Jf1 = value(Jfn_1), Lf1 = value(Lfn_1)
disp('****************************');
disp('****滤波器参数矩阵 模态 2****');
disp('****************************');
Af2 = value(inv(value(Vn_2))*Afn_2), Bf2 = value(inv(value(Vn_2))*Bfn_1), Jf2 = value(Jfn_2), Lf2 = value(Lfn_2)
disp(['运行时间: ', num2str(toc)]); 