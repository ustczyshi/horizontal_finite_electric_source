%% 电场计算和磁场计算的中间项,
% 中间项是对应阶跃响应的中间项
% 对应文档中的中间量C_e1~e11;
% 对应文档中的中间量C_e2~e12;
% 对应文档中的中间量C_h1~h11;
% 对应文档中的中间量C_h2~h12;
function [ e11_1_impulse,e11_01_step,e12_1_impulse,e12_01_step,h11_1_impulse,h11_01_step,h12_1_impulse,h12_01_step] = calculate_temp(r1,r2,z,t,G_S,m2,J_1,delta_1)
u0 = 4*pi*1e-7;
% 电场积分项
% 对应文档中的中间量C_e1~e11;
e11_1_impulse = zeros(1,length(t));
e11_01_step = zeros(1,length(t));
% 对应文档中的中间量C_e2~e12;
e12_1_impulse = zeros(1,length(t));
e12_01_step = zeros(1,length(t));
% 磁场积分项
% 对应文档中的中间量C_h1~h11;
h11_1_impulse = zeros(1,length(t));
h11_01_step = zeros(1,length(t));
% 对应文档中的中间量C_h2~h12;
h12_1_impulse = zeros(1,length(t));
h12_01_step = zeros(1,length(t));

%计算lambda，并将lambda和frequency扩展成二维矩阵-----针对垂直方向汉克尔变换
lambda_r1 = (1./r1) .*exp(delta_1); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda，针对R1接地项
lambda_r2 = (1./r2) .*exp(delta_1); % 如何计算lambda，由采样点的横坐标偏移量转换为积分变量lambda,针对R2接地项

for ii=1:length(t)
     freq = (log(2)*1i/(t(ii)*2*pi))*m2;
     %-------------------------------------------------------
     %根据递推公式计算空气-地面的反射系数
      [lambda1_Array,frequency1_Array] = Array_trans(lambda_r1,freq); % 针对R1接地项
      [lambda2_Array,frequency2_Array] = Array_trans(lambda_r2,freq); % 针对R2接地项  
      r_TE1=calculate_r_TE(lambda_r1,freq); % 针对R1接地项
      r_TE2=calculate_r_TE(lambda_r2,freq); % 针对R2接地项  
      
      z1bar_1 = calculate_r_TM_Z1bar(lambda_r1,freq);%  针对R1接地项  
      z1bar_2 = calculate_r_TM_Z1bar(lambda_r2,freq);%  针对R2接地项  
      %%  计算空气中的阻抗率 zo_bar
      w1 = 2.*pi.*frequency1_Array;%矩阵
      w2 = 2.*pi.*frequency2_Array;%矩阵
      
      z0_bar1 = 1i.*w1.*u0;
      z0_bar2 = 1i.*w2.*u0;
      %% ---------------------------------------------------针对R1接地项的中间值
      % 电场接地项的中间值
      f_e11 = exp(lambda1_Array.*(z)).*(2.*z1bar_1-(1+r_TE1).*(z0_bar1)./lambda1_Array); 
      g_e11 = Fast_Hankel(r1,f_e11,J_1);%列向量
      e11_1_impulse(ii) = GS_Trans2(t(ii),g_e11,G_S); %正脉冲响应时域
      e11_01_step(ii) = GS_Trans(t(ii),g_e11,freq,G_S);%正阶跃响应时域
      % 磁场接地项的中间值
      f_h11 = (1+r_TE1).*exp(lambda1_Array.*(z));% r_TM1 = 1;
      g_h11 = Fast_Hankel(r1,f_h11,J_1);%正阶跃响应频域，f_h11的同一行的对应不同lambda（汉克尔变换滤波系数不同偏移量）的值
      h11_1_impulse(ii) = GS_Trans2(t(ii),g_h11,G_S); %正脉冲响应时域
      h11_01_step(ii) = GS_Trans(t(ii),g_h11,freq,G_S);%正阶跃响应时域
      %% ---------------------------------------------------针对R2接地项的中间值
      % 电场接地项的中间值
      f_e12 = exp(lambda2_Array.*(z)).*(2.*z1bar_2-(1+r_TE2).*(z0_bar2)./lambda2_Array); 
      g_e12 = Fast_Hankel(r2,f_e12,J_1);%列向量
      e12_1_impulse(ii) = GS_Trans2(t(ii),g_e12,G_S); %正脉冲响应时域
      e12_01_step(ii) = GS_Trans(t(ii),g_e12,freq,G_S);%正阶跃响应时域
      % 磁场接地项的中间值
      f_h12 = (1+r_TE2).*exp(lambda2_Array.*(z));% r_TM2 = 1;
      g_h12 = Fast_Hankel(r2,f_h12,J_1);%正阶跃响应频域，f_h11的同一行的对应不同lambda（汉克尔变换滤波系数不同偏移量）的值
      h12_1_impulse(ii) = GS_Trans2(t(ii),g_h12,G_S); %正脉冲响应时域
      h12_01_step(ii) = GS_Trans(t(ii),g_h12,freq,G_S);%正阶跃响应时域
end
% temp_e1  = e11_01_step;
% temp_e2  = e12_01_step;
% temp_h1  = h11_01_step;
% temp_h2  = h12_01_step;


end
