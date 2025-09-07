function dy = Diprotic_ODE(t,y,K)
dy = zeros(6,1);    % a column vector
dy(1) = -K(1)*y(1) +K(2)*y(2)*y(5);
dy(2) = K(1)*y(1)-K(2)*y(2)*y(5)-K(3)*y(2)+K(4)*y(5)*y(3);
dy(3) = K(3)*y(2)-K(4)*y(5)*y(3);
dy(4) = -K(5)*y(4);
dy(5) = K(6)-K(7)*y(5)*y(6)+K(3)*y(2)-K(4)*y(5)*y(3)+K(1)*y(1)-K(2)*y(2)*y(5);
dy(6) = K(6)-K(7)*y(5)*y(6)+K(5)*y(4);

