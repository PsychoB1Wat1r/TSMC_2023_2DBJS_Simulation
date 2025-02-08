%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022-08-13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Metal Rolling Process 后向差分 Auto XieLiHua 2006 %%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
%% 轧制参数
rho1 = 1800; % 弹簧刚度
rho2 = 1800; % 金属硬度
dt = 0.8; % 步长
M = 100; % 辊缝调整机构总质量
rho = (rho1 * rho2) / (rho1 + rho2); % 复合钢度
%% 参数
phi_1 = (2 * M) / (rho * dt^2 + M);
phi_2 = -M / (rho * dt^2 + M);
phi_3 = rho*(dt^2 + M / rho) / (rho * dt^2 + M);
phi_4 = -2 * rho * M / (rho1*(rho * dt^2 + M));
phi_5 = rho * M / (rho1*(rho * dt^2 + M));
b = (-rho * dt^2) / (rho2*(rho * dt^2 + M));

A = [phi_3   phi_4   phi_1   phi_2   phi_5;
     0       0       1       0       0;
     phi_3   phi_4   phi_1   phi_2   phi_5;
     0       0       1       0       0;
     0       1       0       0       0];
B = [b 0 b 0 0]';
%% Feedback Controller Gain
K1 = [1300 0 0 0 0];
K2 = [1200 0 0 0 0];
%% 系统参数
A_Bernoulli_1 = A + B*K1;
A_Bernoulli_2 = A + B*K2;
disp('***********************');
disp('********系统参数********');
disp('***********************');
A1 = roundn(A_Bernoulli_1, -4)
A2 = roundn(A_Bernoulli_2, -4)
eig(A1)
eig(A2)
disp(['运行时间: ', num2str(toc)]); 