function dy = Monoprotic_ODE(t,y,K)
dy = zeros(5,1);    % a column vector

dy(1) = -K(1)*y(1) +K(2)*y(2)*y(4);
dy(2) = +K(1)*y(1) -K(2)*y(2)*y(4);
dy(3) = -K(3)*y(3);
dy(4) = K(1)*y(1) -K(2)*y(2)*y(4)+K(4)-K(5)*y(4)*y(5);
dy(5) = K(4)-K(5)*y(4)*y(5)+K(3)*y(3);

