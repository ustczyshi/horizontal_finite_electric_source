%% ......function description
% function name :Vertical_Magnetic_Dipole_Step_Response
% @u0：真空中的磁导率
% @rou：大地电阻率
% @r：发射线圈半径
% @Tob：观测时间

% 正阶跃和负阶跃不是简单的符号相反的关系
function [step_H01,step_H10,impluse_H10] = Central_Circle_Loop_Step_Response(u0,rou,I,r,Tob)
sita = (u0./rou./Tob).^(0.5)./2;
temp = sita.*r;
Hcoeff = I./2./r;
dHdtcoeff = -I*rou./u0./r^.3;
% dHdtcoeff = -I*rou./r^.3;
step_H10 = Hcoeff.*(erf(temp).*(1-3./2./temp.^2)+3*exp(-temp.^2)./pi^0.5./temp);% 负阶跃
step_H10 = [Hcoeff step_H10(1:length(Tob)-1)];
step_H01 = Hcoeff-step_H10;%正阶跃响应
impluse_H10 = dHdtcoeff.*(3.*erf(temp)-2.*temp./pi.^0.5.*(3+2.*temp.^2).*exp(-temp.^2));% 负脉冲
impluse_H10 = [0 impluse_H10(1:length(Tob)-1)];
end