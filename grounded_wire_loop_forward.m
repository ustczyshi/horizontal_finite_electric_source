%--------------------------------------------------------------------------
%仿真条件：地面TEM 半航空TEM
% 地面水平接地长导线源激励，地面或空中观测 B场

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
x = 100; % 收发水平偏移距，沿y轴
y = 100;
L= 1; % 发射线缆长度，沿x轴
I = 1; % 发射电流
%%  半航空收发高度参数
z =0;% 观测点距地面的高度，地面以上为负值
h =0;% 源距地面的高度

%% 采样率和观测时间段设置
fs = 1e5;% 采样率
dt = 1./fs;
t = 1/fs:1/fs:4e-2;% 时间区间
%%

[hz_01,hz_10,hz_1_impulse,hx_01,hx_10,hx_1_impulse,hy_01,hy_10,hy_1_impulse,ex_01,ex_impulse,ey_01,ey_impulse] = Calculate_Horizontal_Finite_Electrical_Source(I,L,h,x,y,z,t);
% save('horizontal_electrical_dipole_impulse_shuzhijie','hx_1_impulse','hy_1_impulse','hz_1_impulse');
%% Ex-Ey 
%
% Ex -Ey step response
figure;
loglog(t,abs(ex_01) ,'r','linewidth',2);
hold on;
loglog(t,abs(ey_01) ,'b','linewidth',2);
hold on;
grid on;
legend('数值解ex\_1\_step','数值解ey\_1\_step');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop  step response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
% Ex -Ey impulse response
figure;
loglog(t,abs(ex_impulse) ,'r','linewidth',2);
hold on;
loglog(t,abs(ey_impulse) ,'b','linewidth',2);
hold on;
grid on;
legend('数值解ex\_impulse','数值解ey\_impulse');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop  impulse response'])
xlabel('Time/(ms)')
ylabel('E/(V/m)');
%
%%
% hx-hy-hz impulse response
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hz_01)./dt),'b','Linewidth',2);
grid on;
legend('数值解hz\_1\_impulse','数值解diff(hz\_01)./dt');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz impulse response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hx_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hx_01)./dt),'b','Linewidth',2);
grid on;
legend('数值解hx\_1\_impulse','数值解diff(hx\_01)./dt');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx impulse response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hy_1_impulse),'r','Linewidth',2);
hold on
loglog(t(1:end-1).*10^3,u0.*abs(diff(hy_01)./dt),'b','Linewidth',2);
grid on;
legend('数值解hy\_1\_impulse','数值解diff(hy\_01)./dt');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By impulse response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%%
% hx-hy-hz step response
figure;
loglog(t(1:end).*10^3,u0.*abs(hz_01),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(hz_10),'b','Linewidth',2);
grid on;
legend('数值解hz\_01','数值解hz\_10');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bz step response'])
xlabel('Time/(ms)')
ylabel('Bz/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hx_01),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(hx_10),'b','Linewidth',2);
grid on;
legend('数值解hx\_01','数值解hx\_10');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop Bx step response'])
xlabel('Time/(ms)')
ylabel('Bx/(T)');

figure;
loglog(t(1:end).*10^3,u0.*abs(hy_01),'r','Linewidth',2);
hold on
loglog(t(1:end).*10^3,u0.*abs(hy_10),'b','Linewidth',2);
grid on;
legend('数值解hy\_01','数值解hy\_10');
title(['source moment' num2str(I) 'm*' num2str(L) 'm position ('  num2str(x) ',' num2str(y) ',' num2str(z) ')wire loop By step response'])
xlabel('Time/(ms)')
ylabel('By/(T)');
%% save data

