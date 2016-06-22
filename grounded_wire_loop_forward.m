%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，地面或空中观测 B场
% 1. 均匀半空间电阻率设定100欧姆米,先观测不同高度响应的变化，主要是观察Hz,Hx,Hy,高度设定0,30,60,90米
% 
%--------------------------------------------------------------------------
%%
format long;
clear all;clc;close all;
%%
u0 = 4*pi*1e-7;
load parameters.txt;
sigma1 = parameters(1,2);%第一层的电导率
rou = 1./sigma1;
%% 发射机参数地面
x = 0; % 收发水平偏移距，沿y轴
y = 100;
L= 1000; % 发射线缆长度，沿x轴
I = 1; % 发射电流
m = I*L;
%%  半航空收发高度参数
z =0;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度

%% 采样率和观测时间段设置
fs = 1e5;% 采样率
dt = 1./fs;
t = 1/fs:1/fs:4e-2;% 时间区间
%%

%% 解析解
% 有限长接地导线中垂线地面上的解析解
% [ hz_impulse_jiexijie] = ground_finite_wire_source_jiexi(u0,rou,I,L,y,t);
% save(['ground_finite_wire_source_jiexi_y' num2str(y) '.mat'],'hz_impulse_jiexijie');
% 水平电偶源
% [step_ex01,step_ey01,impluse_hx,impulse_hy,impulse_hz] = Horizontal_Electrical_Dipole_Response(u0,rou,m,x,y,z,t);

[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] =...
            Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t);
        % [ hz_impulse_jiexijie] = ground_finite_wire_source_jiexi(u0,rou,I,L,y,t);
%% Ex-Ey 
%{
% 当接地长导线尺寸源小于偏移距是，长导线瞬变电磁正演与水平电偶源响应的对比
% Ex  step response
figure;
loglog(t,abs(ex_01) ,'r','linewidth',2);
hold on;
loglog(t,abs(step_ex01) ,'b','linewidth',1);
hold on;
grid on;
legend('数值解ex\_1\_step','解析解step\_ex01');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop  Ex step response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
% Ey impulse response
figure;
loglog(t,abs(ey_01) ,'r','linewidth',2);
hold on;
loglog(t,abs(step_ey01) ,'b','linewidth',1);
hold on;
grid on;
legend('数值解ey_01','解析解step_ey01');
title(['source moment' num2str(I) 'Ax' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Ey step response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
%
% Hz impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impulse_hz),'b','Linewidth',1);
grid on;
legend('数值解hz\_1\_impulse','解析解impulse\_hz');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');
% Hx impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hx_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impluse_hx),'b','Linewidth',1);
grid on;
legend('数值解hx\_1\_impulse','解析解impluse\_hx');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx step response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');
% Hy impulse
figure;
loglog(t(1:end).*10^3,u0.*abs(hy_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(impulse_hy),'b','Linewidth',1);
grid on;
legend('数值解hy\_1\_impulse','解析解impulse\_hy');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By step response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%}
%%
% hz impulse response
% figure;
% loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
% hold on
% loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./dt),'b','Linewidth',2);
% grid on;
% legend('数值解hz\_1\_impulse','数值解diff(hz\_01)./dt');
% title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
% xlabel('Time/(ms)')
% ylabel('Bz/(T)');
%{
%% 接地长导线的hz数值解和解析解对比
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(hz_impulse_jiexijie),'b','Linewidth',2);
grid on;
legend('数值解hz\_1\_impulse','解析解hz\_impulse\_jiexijie');
title(['source moment' num2str(I) 'A*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');
%%
% 验证接地长导线的Bz响应
figure;
plot(t(1:end).*1e3,abs(hz_impulse_jiexijie)./abs(hz_1_impulse));
title('Bz数值解和解析解的误差')
xlabel('Time/(ms)')
ylabel('error/(T)');

figure;
plot(t(1:end).*1e3,(abs(hz_impulse_jiexijie)./abs(hz_1_impulse)-1).*100);
title('Bz数值解和解析解的误差')
xlabel('Time/(ms)')
ylabel('error/(%)');
%}



