close all;
L=100E-9;
zeta_max=5E-9;
L=500;
zeta_max=20;
zeta_max=L/8;
X=linspace(0,L,500);
w=zeta_max/8;
zeta_max/w
%Y=zeta_max*(X/L-sin(2*acos(-1.)*X/L)/(2*acos(-1.)));
Y=1*exp(-(X-zeta_max).^2/(2*w^2));
%Y=1*exp(-(X-L/5).^2/(2*(L/(20))^2));
%Y=1*exp(-(X-L/5).^2*200/L.^2)
%Y=1*exp((200/L).*(1/5-X./L))
plot(X,Y,2*zeta_max*[1 1],max(Y)*[0 1])