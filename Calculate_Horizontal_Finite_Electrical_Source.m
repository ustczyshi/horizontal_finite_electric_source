% %% 计算水平电偶极源各积分项的结果
% function [hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse] = Calculate_Horizontal_Electrical_Dipole(I,L,h,x,y,z,t)
%sigma1 :第一层的电导率
% 调用子函数
%{
[ temp_e1,temp_e2,temp_h1,temp_h2] = calculate_temp(r1,r2,z,t,G_S,m2,J_1,delta_1)
%}

function [hz_01,hz_10,hz_impulse,hx_01,hx_10,hx_impulse,hy_01,hy_10,hy_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t)
u0 = 4*pi*1e-7;
r = (x^2+y^2).^(0.5);
%% 计算接地项的参数
R1 = ((x+L/2).^2+y.^2).^0.5; 
R2 = ((x-L/2).^2+y.^2).^0.5;

%% 磁场积分项时域表达
hz1_1_impulse = zeros(1,length(t));% 正脉冲响应时域
hz1_01_step = zeros(1,length(t));% 正阶跃响应时域
hz1_10_step = zeros(1,length(t));% 正阶跃一次场

h1_1_impulse = zeros(1,length(t));% 正脉冲响应时域
h1_01_step = zeros(1,length(t));% 正阶跃响应时域
% h1_01_step1 = zeros(1,length(t));% 正阶跃一次场

h01_1_impulse = zeros(1,length(t));% 正脉冲响应时域
h01_01_step = zeros(1,length(t));% 正阶跃响应时域
% h01_01_step1 = zeros(1,length(t));% 正阶跃一次场

h00_1_impulse = zeros(1,length(t));% 正脉冲响应时域
h00_01_step = zeros(1,length(t));% 正阶跃响应时域
% h00_01_step1 = zeros(1,length(t));% 正阶跃一次场
hy_01_step1= zeros(1,length(t));% 正阶跃一次场
% 电场积分项
e11_1_impulse = zeros(1,length(t));
e11_01_step = zeros(1,length(t));
e20_1_impulse = zeros(1,length(t));
e20_01_step = zeros(1,length(t));
e0_1_impulse = zeros(1,length(t));
e0_01_step = zeros(1,length(t));
%% 
G_S=load ('G_S.txt')';% G_S行向量
m2 = 1:length(G_S);
%--------------------------------------------------------------------------
%第一步：读取已经存储的滤波器系数,表示为行向量；
%%  计算e0,h00,h01
load J0_Gupt.txt;       
J_zero = J0_Gupt( :, 3)'; % 快速汉克尔变换滤波系数
delta = J0_Gupt( :, 2)'; %  采样点的横坐标偏移量
%% e1 h1 方向计算
 load J1_Gupt.txt;       
J_1 = J1_Gupt( :, 3)'; % 快速汉克尔变换滤波系数
delta_1 = J1_Gupt( :, 2)'; %  采样点的横坐标偏移量
% 计算中间值
[ e11_impulse,e11_step,e12_impulse,e12_step,h11_impulse,h11_step,h12_impulse,h12_step] = calculate_temp(R1,R2,z,t,G_S,m2,J_1,delta_1);
% 计算Ey-Hx
ey_01 =real(- I.*y./(4.*pi).*(1./R2.*e12_step-1./R1.*e11_step));
ey_impulse = real(-I.*y./(4.*pi).*(1./R2.*e12_impulse-1./R1.*e11_impulse));% 正阶跃响应
% ey_10 = -ey_01;
hx_01 = I.*y./(4.*pi).*(1./R2.*h12_step-1./R1.*h11_step);% + step response,only the second response
 hx_10 = -hx_01;% -step response,only the second response
hx_impulse = I.*y./(4.*pi).*(1./R2.*h12_impulse-1./R1.*h11_impulse);% + step response,only the second response 
% 计算Ex-Hy的与接地项有关的部分




%% -----------------------------------------计算沿长导线的积分项---------------------------------------------------------------
%计算lambda，并将lambda和frequency扩展成二维矩阵-----针对垂直方向汉克尔变换
lambda_0 = (1./r) .*exp(delta); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda，针对J0
lambda_1 = (1./r) .*exp(delta_1); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda,针对J1汉克尔变换
 for ii=1:length(t)
        freq = (log(2)*1i/(t(ii)*2*pi))*m2;
        %--------------------------------------------------------------------------
        %根据递推公式计算空气-地面的反射系数
        [lambda0_Array,frequency0_Array] = Array_trans(lambda_0,freq); % 针对针对J0汉克尔变换
        [lambda1_Array,frequency1_Array] = Array_trans(lambda_1,freq); % 针对针对J1汉克尔变换   
        r_TE0=calculate_r_TE(lambda_0,freq); % 针对J0
        r_TE1=calculate_r_TE(lambda_1,freq); % 针对J1
        z1bar_1 = calculate_r_TM_Z1bar(lambda_1,freq);% 针对J1
        z1bar_0 = calculate_r_TM_Z1bar(lambda_0,freq);% 针对J0
        %垂直磁偶源z轴磁场的核函数,h是源位置,z为观测位置
        % 垂直磁偶极源z轴磁场的核函数 u_0 = lambda,z=h=0,地面发射地面接收;
        %% -----------------------------------半航空tem 电场分量积分项 h=0----------------------------------------
        w0 = 1i.*2.*pi.*frequency0_Array;%矩阵
        w1 = 1i.*2.*pi.*frequency1_Array;%矩阵
         w =  1i.*2.*pi.*freq';
        %% 针对Ex
        f_e0 = (1+r_TE0).*exp(lambda0_Array.*(z));
        g_e0 = Fast_Hankel(r,f_e0,J_zero);
        e0_1_impulse(ii) = GS_Trans2(t(ii),w.*g_e0,G_S); %正脉冲响应时域
        e0_01_step(ii) = GS_Trans(t(ii),w.*g_e0,freq,G_S);%正阶跃响应时域
        
        %% -----------------------------------半航空tem 磁场分量积分项 h=0----------------------------------------
        %% 针对Hz
            f_hz1 = (1+r_TE1).*exp(lambda1_Array.*(z)).*lambda1_Array;% r_TM1 = 1;
            primary_hz1 = exp(lambda1_Array.*(z)).*lambda1_Array;%% 对应直流成分，此时的核函数中的r_TE=0  
            g_hz1 = Fast_Hankel(r,f_hz1,J_1);%正阶跃响应频域
            g_hz1_zeros = Fast_Hankel(r,primary_hz1,J_1);%正阶跃响应零频响应
%         h1_0_impulse(ii) = -GS_Trans2(t(ii),g_h1,G_S);%负脉冲响应时域
%         h1_10(ii) = +GS_Trans(t(ii), g_h1_zeros,freq,G_S)-GS_Trans(t(ii),g_h1,freq,G_S);%负阶跃响应时域
            hz1_10_step(ii) = +GS_Trans(t(ii), g_hz1_zeros,freq,G_S)-GS_Trans(t(ii),g_hz1,freq,G_S);%负阶跃响应时域
            hz1_1_impulse(ii) = GS_Trans2(t(ii),g_hz1,G_S); %正脉冲响应时域
            hz1_01_step(ii) = GS_Trans(t(ii),g_hz1,freq,G_S);%正阶跃响应时域

       %% 针对 Hy
       f_h00 = (1-r_TE0).*exp(lambda0_Array.*(z)).*lambda0_Array;
       primary_h00 =-exp(lambda0_Array.*(z)).*lambda0_Array;%% 对应直流成分，此时的核函数中的r_TE=0
       g_h00 = Fast_Hankel(r,f_h00,J_zero);%正阶跃响应频域
       g_h00_zeros = Fast_Hankel(r,primary_h00,J_zero);%正阶跃响应零频响应
%        h00_0_impulse(ii) = -GS_Trans2(t(ii),g_h00,G_S);%负脉冲响应时域
       hy_01_step1(ii) = +GS_Trans(t(ii), g_h00_zeros,freq,G_S);% 一次场
%        h00_10(ii) = +GS_Trans(t(ii), g_h00_zeros,freq,G_S)-GS_Trans(t(ii),g_h00,freq,G_S);%负阶跃响应时域
       h00_1_impulse(ii) = GS_Trans2(t(ii),g_h00,G_S); %正脉冲响应时域
       h00_01_step(ii) = GS_Trans(t(ii),g_h00,freq,G_S);%正阶跃响应时域
 end
    %电场分量
        ex_01 =  real(-I.*L./(4.*pi).*((y.^2-x.^2)./r.^3.*e11_01_step+(x.^2)./r.^2.*e20_01_step)-I.*L.*u0./(4.*pi).*e0_01_step);
        ex_impulse =real(- I.*L./(4.*pi).*((y.^2-x.^2)./r.^3.*e11_1_impulse+(x.^2)./r.^2.*e20_1_impulse)-I.*L.*u0./(4.*pi).*e0_1_impulse);  
    % 磁场分量
        hz_01 = I.*L.*y./(4.*pi.*r).*hz1_01_step;% + step response
        hz_10 = I.*L.*y./(4.*pi.*r).*hz1_10_step;%- step response
        
        hy_01 = -I.*L./(4.*pi).*((y.^2-x.^2)./r.^3.*h1_01_step+x.^2./r.^2.*h01_01_step+ h00_01_step);% + step response
         hy_10 = I.*L./(4.*pi).*hy_01_step1 - hy_01;
%         hy_10 = hy_01(end) - hy_01;
      
        hz_impulse = I.*L.*y./(4.*pi.*r).* hz1_1_impulse;
        hy_impulse = -I.*L./(4.*pi).*((y.^2-x.^2)./r.^3.*h1_1_impulse+x.^2./r.^2.*h01_1_impulse+h00_1_impulse);

        load parameters.txt;% 读取参数
       save(['semiatem_horizontal_electrical_dipole_response_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y) '.mat'],...
        'hz_01','hz_10','hz_impulse','hx_01','hx_10','hx_impulse','hy_01','hy_10','hy_impulse','ex_01','ex_impulse','ey_01','ey_impulse','parameters');
end