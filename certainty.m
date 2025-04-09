%% 微电网两阶段鲁棒优化经济调度方法

clc;
clear;
yalmip('clear')


%% 参数
% 微型燃气轮机
    P_Gmin = 80;
    P_Gmax = 800;
    a = 0.67;

% 储能
    P_Smax = 500;
    E_Smax = 1800;
    E_Smin = 400;
    E_S0 = 1000;
    K_S = 0.38;
    yita = 0.95;
% 需求响应负荷
    K_DR = 0.32;
    D_DR = 2940;
    DR_max = 200;
    DR_min = 50;
    P_DRstar = ...
    [   8.2e+001
        7.1e+001
        6.1e+001
        5.1e+001
        7.1e+001
        7.3e+001
        9.1e+001
        1.02e+002
        1.22e+002
        1.54e+002
        1.70e+002
        2.00e+002
        1.39e+002
        1.03e+002
        1.01e+002
        1.22e+002
        1.40e+002
        1.50e+002
        1.90e+002
        2.00e+002
        2.00e+002
        1.91e+002
        1.01e+002
        8.2e+001];
% 配电网交互功率
    P_Mmax = 1500;

%% 不确定集 [u_pv, u_L]
u_pv_est = [
    0
    0
    0
    0
    0
    4.17074232453333e+000
    5.94883875758584e+001
    2.48146877320564e+002
    5.36796660947943e+002
    8.81018103602796e+002
    1.09189870176427e+003
    1.18055770363352e+003
    1.13031213991243e+003
    9.30100120349269e+002
    8.79862750620951e+002
    8.24069853788441e+002
    5.90516477607354e+002
    3.95835403170051e+002
    1.28924282385476e+002
    0
    0
    0
    0
    0];

u_L_est = [
    3.47699739155637e+002
    3.24362863491450e+002
    3.08909284127137e+002
    2.89508260071491e+002
    3.37179016520143e+002
    4.67769297652401e+002
    5.39194280745822e+002
    5.55273886581007e+002
    6.06834122307023e+002
    7.41418220461791e+002
    8.44399574920297e+002
    9.11819147908415e+002
    7.38490967056323e+002
    6.95405274852671e+002
    6.79951695488359e+002
    7.27634045019805e+002
    7.71392136025505e+002
    7.91442372717612e+002
    8.43048980774804e+002
    8.90800888803014e+002
    8.35838083276978e+002
    5.75735677712298e+002
    4.61561201816250e+002
    4.02639358516085e+002];

%% 设置价值向量
lambda_t = [
    0.45
    0.45
    0.45
    0.45
    0.45
    0.45
    0.45
    0.9
    1.35
    1.35
    1.35
    0.9
    0.9
    0.9
    0.9
    0.9
    0.9
    0.9
    1.35
    1.35
    1.35
    1.35
    1.35
    0.45];  % 日前交易电价

%% 约束矩阵

E24 = eye(24,24);  % 24*24 单位矩阵
Z24 = zeros(24,24);    % 24*24 零矩阵
L24 = tril(ones(24,24));    % 24*24 下三角形为1阵

% 目标函数: c'*y
% 约束(1),(3),(11),(17)
    c = [
               a*ones(24,1);
        K_S*yita*ones(24,1);
        K_S/yita*ones(24,1);
                zeros(24,1);
            K_DR*ones(48,1);
                   lambda_t;
                  -lambda_t;
                zeros(48,1);
        ];

% D*y >= d       
% 约束(2),(7),(9),(13)
    D = [
        E24,        Z24,           Z24,       Z24,    Z24,  Z24;
       -E24,        Z24,           Z24,       Z24,    Z24,  Z24;
        Z24,  yita.*L24,  -1/yita.*L24,       Z24,    Z24,  Z24; 
        Z24, -yita.*L24,   1/yita.*L24,       Z24,    Z24,  Z24;
        Z24,        Z24,           Z24,       E24,    Z24,  Z24;
        Z24,        Z24,           Z24,      -E24,    Z24,  Z24;
        Z24,        Z24,           Z24,       Z24,    E24,  Z24;
        Z24,        Z24,           Z24,       Z24,    Z24,  E24;
        ];

    D = [D, zeros(24*8, 24*4)];

    d = [
            P_Gmin*ones(24,1);
            -P_Gmax*ones(24,1);
            (E_Smin - E_S0)*ones(24,1);
            (E_S0 - E_Smax)*ones(24,1);
            DR_min*ones(24,1);
            -DR_max*ones(24,1);
            zeros(24,1);
            zeros(24,1);
        ];


% K y == g 
% 约束(6),(8),(12),(14)
    K = [
         zeros(1,24),   yita*ones(1,24),   -1/yita*ones(1,24), zeros(1,24*7);
         zeros(1,24*3), ones(1,24),     zeros(1,24*6);
         Z24,               Z24,                  Z24,          E24,   E24, -E24,  Z24,  Z24,  Z24,  Z24;
        -E24,               E24,                 -E24,          E24,   Z24,  Z24, -E24,  E24, -E24,  E24;
        ];

    g = [
        0;
        D_DR;
        P_DRstar;
        zeros(24,1)
        ];
 
% Fx + Gy >= h  
% 约束(4),(5),(15),(16)
    F = [
         P_Smax*E24,            Z24;
        -P_Smax*E24,            Z24;
                Z24,     P_Mmax*E24;
                Z24,    -P_Mmax*E24;
        ];

    G = [
        Z24,  Z24, -E24, Z24, Z24, Z24,  Z24,  Z24, Z24, Z24;
        Z24, -E24,  Z24, Z24, Z24, Z24,  Z24,  Z24, Z24, Z24;
        Z24,  Z24,  Z24, Z24, Z24, Z24, -E24,  Z24, Z24, Z24;
        Z24,  Z24,  Z24, Z24, Z24, Z24,  Z24, -E24, Z24, Z24;
        ];

    h = [
               zeros(24,1);
        -P_Smax*ones(24,1);
               zeros(24,1);
        -P_Mmax*ones(24,1);
        ];


% I_u y == u_hat 
    I_u = [
        zeros(24,192), E24, Z24;
        zeros(24,192), Z24, E24;
        ];


%%  变量
% 第一阶段决策变量 x = [U_S; U_M]
U_S = binvar(24,1);
U_M = binvar(24,1);
x = [U_S; U_M];
%%% 第二阶段决策变量 y = [P_G; P_Sch; P_Sdis; P_DR; P_DR1; P_DR2; P_Mbuy; P_Msell; P_pv; P_L];
P_G = sdpvar(24,1);
P_Sch = sdpvar(24,1);
P_Sdis = sdpvar(24,1);
P_DR = sdpvar(24,1);
P_DR1 = sdpvar(24,1);
P_DR2 = sdpvar(24,1);
P_Mbuy = sdpvar(24,1);
P_Msell = sdpvar(24,1);
P_pv = sdpvar(24,1);
P_L = sdpvar(24,1);
y = [P_G; P_Sch; P_Sdis; P_DR; P_DR1; P_DR2; P_Mbuy; P_Msell; P_pv; P_L];

%% 循环变量
% 迭代次数
itr = 10;

% u_k存储
u_pv_k = zeros(24,itr);
u_L_k = zeros(24,itr);
y_k = sdpvar(240,itr,'full');
alpha = sdpvar(1);

% 不确定变量
u_pv = sdpvar(24,1);
u_L = sdpvar(24,1);
u = [u_pv; u_L];
u_hat = [u_pv_est; u_L_est];
du = [-u_pv_est*0.15; u_L_est*0.1];

% 对偶变量
gamma = sdpvar(size(D,1),1);
lamda = sdpvar(size(K,1),1);
nu = sdpvar(size(G,1),1);
pai = sdpvar(size(I_u,1),1);
miu = sdpvar(size(y,1),1);

% 不确定 01变量
B_pv = binvar(24,1);
B_L = binvar(24,1);
B = [B_pv;B_L];



%% loop
% 上下界
UB = inf;
LB = -inf;
lb = []; ub = [];

% 大 M 法用到的 M
M = 1e5;

% 不确定性调节参数
Gamma_pv = 6;
Gamma_L = 12;

% 求解器设置
opt = sdpsettings('verbose', 0 ,'solver','gurobi');
% ops.gurobi.MIPGap = 1e-3; 
% ops.gurobi.TimeLimit = 5;      % 设置求解时间为10s

Cons = [];
Obj = c'*y;

Cons = [
    D*y >= d;
    K*y == g;
    F*x + G*y >= h;
    I_u*y == u_hat;
    y >= 0;
    ];
optimize(Cons, Obj, opt);

obj=value(Obj);

disp("------------------------")
disp(['确定性模型的目标函数值为： ',num2str(obj)])
disp("------------------------")




