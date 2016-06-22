%% 比较相同偏移距下，不同接收高度响应的衰减情况
clc;close all;
%
u0 = 4*pi*1e-7;
fs = 1e5;
dt = 1/fs;
x = 0;
L = 1000;
h = 0;
x = 0;
rou =100;
% y = [100,500,1000,2000];
y = 1000;
z = [0,-30,-60,-90,-150];
% N = 4e-2*fs;
% t = (1:N).*dt;
t = logspace(-5,0,50);% 时间区间
N = length(t);
%% dataname


hz_step_y100 = zeros(length(z),N);
hz_impulse_y100 = zeros(length(z)+1,N);

hy_step_y100 = zeros(length(z),N);
hy_impulse_y100 = zeros(length(z),N);

ex_step_y100 = zeros(length(z),N);
ex_impulse_y100 = zeros(length(z),N);
for k = 1:length(y)
    dataname0 = ['ground_finite_wire_source_jiexi_rou' num2str(rou) '_L' num2str(L) '_y' num2str(y(k)) '.mat'];
    load(dataname0);
    hz_impulse_y100(length(z)+1,:) = -hz_impulse_jiexijie;
    for kk = 1:length(z)
        dataname =['SemiAtem_Horizontal_Finite_Electrical_Source_L' num2str(L) '_h' num2str(h) '_z' num2str(z(kk)) '_x' num2str(x) '_y' num2str(y(k)) '.mat'];
        load(dataname);
        hz_step_y100(kk,:) = hz_10;
        hz_impulse_y100(kk,:) = hz_impulse;
        hy_step_y100(kk,:) = hy_10;
        hy_impulse_y100(kk,:) = hy_impulse;
        ex_step_y100(kk,:) = ex_01;
        ex_impulse_y100(kk,:) = ex_impulse;
    end
end
%% display;
%% Hz
figure;
loglog(t,hz_step_y100,'linewidth',2);
% title('负阶跃响应hz随接收高度的变化');
title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hz/(A/m)');
% percentage difference
figure;
semilogx(t,(1-(hz_step_y100(1:length(z),:))./repmat(hz_step_y100(1,:),length(z),1)).*100,'linewidth',2);
title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hz percentage difference/(%)');
%% 脉冲响应
figure;
loglog(t,(hz_impulse_y100),'linewidth',2);
% title('脉冲响应hz随接收高度的变化');
title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm'],[ '接收高度' num2str(z(1)) 'm的解析解']);
grid on;
xlabel('Time/s');
ylabel('dHz/dt/(A/m)');
% percentage difference
figure;
semilogx(t,(1-(hz_impulse_y100(1:length(z),:))./repmat(hz_impulse_y100(1,:),length(z),1)).*100,'linewidth',2);
title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm'],[ '接收高度' num2str(z(1)) 'm的解析解']);
grid on;
xlabel('Time/s');
ylabel('dHz/dt percentage difference/(%)');
%%
% Hy
figure;
loglog(t,abs(hy_step_y100),'linewidth',2);
title('负阶跃响应hy随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');
figure;
loglog(t,abs(hy_impulse_y100),'linewidth',2);
title('脉冲响应hy随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');
%% Ex
figure;
loglog(t,abs(ex_step_y100),'linewidth',2);
title('正阶跃响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
figure;
loglog(t,abs(ex_impulse_y100),'linewidth',2);
title('脉冲响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
