%% input variable 
function y=original_krig(X)
x1=X(1);
x2=X(2);

y1=(x2-(5.1/(4*pi^2))*x1.^2+5/pi*x1-6).^2+10*(1-1/(8*pi))*cos(x1)+10;
y2=5*cos(0.25*pi*x1).*x2+10*x1;
y3=5*sin(2*x1-0.5*pi).*x2-10*x2.^2;

y=[y1 y2 y3];
%%


