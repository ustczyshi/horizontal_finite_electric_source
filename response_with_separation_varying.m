%% 比较相同偏移距下，不同接收高度响应的衰减情况
clc;close all;
%
u0 = 4*pi*1e-7;
fs = 1e5;
dt = 1/fs;
x = 0;
L = 1;
h = 0;
rou =4;
y = [100,200,400,800,1600];
z =0;
% N = 4e-2*fs;
% t = (1:N).*dt;
t = logspace(-4,1,100);% 时间区间
N = length(t);
%% dataname


hz_step_mat = zeros(length(y),N);
hz_impulse_mat = zeros(length(y),N);
hz_impulse_jiexi = zeros(length(y),N);
hy_step_mat = zeros(length(y),N);
hy_impulse_mat = zeros(length(y),N);

ex_step = zeros(length(y),N);
ex_impulse_mat = zeros(length(y),N);
for kk = 1:length(y)
    dataname0 = ['ground_finite_wire_source_jiexi_rou' num2str(rou) '_L' num2str(L) '_y' num2str(y(kk)) '.mat'];
    load(dataname0);
    hz_impulse_jiexi(kk,:) = -hz_impulse_jiexijie;
%     for kk = 1:length(z)
%         dataname =['SemiAtem_Horizontal_Finite_Electrical_Source_L' num2str(L) '_h' num2str(h) '_z' num2str(z(kk)) '_x' num2str(x) '_y' num2str(y(k)) '.mat'];
        dataname = ['SemiAtem_Horizontal_Finite_Electrical_Source_separation_varying_L'...
            num2str(L) '_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y(kk)) '.mat'];
%         dataname = ['SemiAtem_Horizontal_Electrical_Dipole_separation_varying_L'...
%             num2str(L) '_h' num2str(h) '_z' num2str(z) '_x' num2str(x) '_y' num2str(y(kk)) '.mat'];
        load(dataname);
        hz_step_mat(kk,:) = hz_10;
        hz_impulse_mat(kk,:) = hz_impulse;
        hy_step_mat(kk,:) = hy_10;
        hy_impulse_mat(kk,:) = hy_impulse;
        ex_step(kk,:) = ex_01;
        ex_impulse_mat(kk,:) = ex_impulse;
%     end
end
%% display;
% Hz
figure;
loglog(t,hz_step_mat,'Linewidth',2);
title('不同偏移距下负阶跃响应Hz');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hz/(A/m)');
%%  percentage difference
figure;
loglog(t(1:end),u0.*hz_impulse_mat,'Linewidth',1);%./repmat(diff(t),length(y),1)
hold on;
loglog(t(1:end),u0.*hz_impulse_jiexi,'--','Linewidth',1);
hold on;
title(['均匀半空间 rou=' num2str(rou)  '不同偏移距下感应电压Uz']);
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Uz/(V/m)');

%%
figure;
% loglog(t(1:end),u0.*hz_impulse,'Linewidth',2);%./repmat(diff(t),length(y),1)
% hold on;
loglog(t(1:end),u0.*hz_impulse_mat,'Linewidth',2);
hold on;
title('不同偏移距下感应电压Uz');
% title(['Grounded Finite wire L=' num2str(L) ' offset=' num2str(y) ' rou =' num2str(rou)]);
legend([  num2str(y(1)) 'm'],[ num2str(y(2)) 'm'],...
    [ num2str(y(3)) 'm'],[ num2str(y(4)) 'm'],[ num2str(y(5)) 'm'],[ num2str(y(1)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Uz/(V/m)');

%{
%%
% Hy
figure;
loglog(t,abs(hy_step_y100));
title('负阶跃响应hy随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');
figure;
loglog(t,abs(hy_impulse_y100));
title('脉冲响应hy随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('Hy/(A/m)');
%% Ex
figure;
loglog(t,abs(ex_step_y100));
title('正阶跃响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
figure;
loglog(t,abs(ex_impulse_y100));
title('脉冲响应ex随接收高度的变化');
legend([ '接收高度' num2str(z(1)) 'm'],[ '接收高度' num2str(z(2)) 'm'],...
    [ '接收高度' num2str(z(3)) 'm'],[ '接收高度' num2str(z(4)) 'm'],[ '接收高度' num2str(z(5)) 'm']);
grid on;
xlabel('Time/s');
ylabel('ex/(V/m)');
%}