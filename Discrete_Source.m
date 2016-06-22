%%    积分离散化为求和
% 这里主要是根据收发偏移距和发射源尺寸，合理设计离散化的最小尺度
% 源的位置在（-L/2,L/2）之间；
% 观测点在（x,y,z）;
% 观测点距源的水平垂直距离为|y|
% 观测点距源离散化后各电偶极的的水平距离为
% min(|y|,((x-L/2)^2+y^2)^(0.5),((x+L/2)^2+y^2)^(0.5))
% |y| 对应观测点x在(-L/2,L/2)之间
% ((x-L/2)^2+y^2)^(0.5),((x+L/2)^2+y^2)^(0.5) 对应观测点的位置在|x|>L/2

%% 离散环的标准
% 当收发距离超过源长度的10倍是可以认为源为电偶极源
% 这里取微元dx尺度为观测点距发射长导线最短距离的1/20为标准
% dx = dmin /20;
% 离散后形成的电偶源的个数N = L/dx;
% r 返回个电偶极源中心距观测点的距离
% 细化因子：1/2~2
function [N, dx,x_center,dmin,r] = Discrete_Source(x,y,z,L,fine)
    r1 = ((x+L/2)^2+y^2)^(0.5);
    r2 = ((x-L/2)^2+y^2)^(0.5);
    %% 计算观测观测点距长导线源的最短距离
    if x >L/2
        dmin = r2;
    elseif x < -L/2
        dmin = r1;
    else
        dmin = abs(y);
    end
    if fine>5 || fine <0.5
        error('细化因子输入错误');
        return;
    end
    zoomout = ceil(50*fine); 
if dmin >= zoomout*L
         r =  (x^2+y^2).^(0.5);
         dx = L;
         N = 1;
         x_center = 0;
         display('可以认为激励源为水平电偶极源，无需分割成小微元');
         return;
     else
        % if the distance from the observation point to the source < 20 times the long of the source
    dx = dmin ./ zoomout;
end
%     if dx > 1
%         dx = 1;
    N =floor(ceil(L./dx)./2).*2+1 ;% the moment cell numbers
    dx = L/N; % real scale of the moment cell
    % 变为奇数，为使中间的电偶源中点落在源点处，也是为了对称考虑
%     x_p = ones(1,(N-1)./2);% locate + x axis
    x_p = dx.*(1:(N-1)/2); % the center position of the moment cell located + x axis direction   
    x_n = -x_p;% the center position of the moment cell located - x axis direction
    x_center = [x_n, 0 , x_p];% the center position of the moment cell
    r = ((x-x_center).^2+y.^2).^0.5;
    %% 绘制结果图看分解效果
    figure;
    x_axis = (-L/2:dx:L/2);
    y_axis = zeros(1,length(x_axis));
    z_axis = zeros(1,length(x_axis));
    y_center = zeros(1,length(x_center));
    z_center = zeros(1,length(x_center));
    plot3(x_axis,y_axis,z_axis,'r','Linewidth',1);% 长导线源的范围
    hold on;
    plot3(x_center,y_center,z_center,'bx','linewidth',1);% 长导线源分成的电偶极源的中心位置
    hold on;
    plot3(x,y,z,'bo','linewidth',1);
    axis([-L L -2*L 2*L -10 10 ]);
    grid on;
    title('simulate model');
    xlabel('x-source scale/m');
    ylabel('y-offset/m');
    zlabel(' z-hight of the receiver/m');
%     end
    
end