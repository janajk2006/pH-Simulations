function out = Monoprotic_DAE_func(t,y,K)

out = [ -K(1)*y(1)+K(2)*y(2)*y(3);
	K(1)*y(1)-K(2)*y(2)*y(3);
	K(3)-K(4)*y(3)*y(4)+K(1)*y(1)-K(2)*y(2)*y(3);
	K(6)-y(4)+y(3)-y(2)];
