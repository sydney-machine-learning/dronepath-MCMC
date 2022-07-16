syms x y z v1 v2 v3 
ntime=1000000;stepsize = 0.005; 
X = [10 40 5 0.1 0.1 0.1]; tau=[45 5 45];

o=[25 25 20]; rad=3;
delta = 17.1972;
delta1=delta;delta2=delta;delta3=delta; beta=0.2568; 

W=0.5*((x-o(1))^2+(y-o(2))^2+(z-o(3))^2-(rad)^2);
V=0.5*((x-tau(1))^2+(y-tau(2))^2+(z-tau(3))^2+v1^2+v2^2+v3^2); 
F=0.5*((x-tau(1))^2+(y-tau(2))^2+(z-tau(3))^2);

L=V+F*beta/W;
sigma1 = -(delta1*v1+diff(L,x));
sigma2 = -(delta2*v2+diff(L,y));
sigma3 = -(delta3*v3+diff(L,z)); 
 
% Draw trajectory
h = plot3(X(1),X(2),X(3),'-');hold on

% Draw target
 ht = plot3(tau(1),tau(2),tau(3),' . ');set(ht,'MarkerSize',10);
hold on 

% Draw Obsatcle
[X1,Y1,Z1] =  sphere(20); 
 hobs = surf(o(1)+X1*rad,o(2)+Y1*rad,o(3)+Z1*rad);set(hobs,'MarkerSize',1, 'FaceColor','r');

% Button to early stop
hstop = uicontrol('Style','pushbutton','String','Stop', 'Position',[30 0 80 10],'callback','earlystop = 1;'); 
earlystop = 0;

axis([0 50 0 50 0 50]) 
grid on 

delete('XYZ.txt')

 for t=1:ntime
     drawnow 
var={x,y,z,v1,v2,v3};
varval={X(1),X(2),X(3),X(4),X(5),X(6)};

% ODEs
FX(1)=X(4);
FX(2)=X(5);
FX(3)=X(6);
FX(4)=subs(sigma1,var,varval); 
FX(5)=subs(sigma2,var,varval);
FX(6)=subs(sigma3,var,varval);
w = (subs(W,var,varval));
stepsize = 0.03*min(1 , 0.02*exp(w/25)+0.5*exp(0.15*(w-16))); 

% RK4 to numerically solve the ODE
 for I = 1:length(X)
kX1(I) = stepsize * FX(I); X(I) = X(I) + kX1 (I) / 2;
kX2(I) = stepsize * FX(I); X(I) = X(I) + kX2 (I) / 2;
kX3(I) = stepsize * FX(I); X(I) = X(I) + kX3 (I) / 2;
X(I) = X(I) + (kX1(I) + 2 * (kX2(I) + kX3(I)) + stepsize * FX(I)) / 6;
 end
 
 if earlystop
    t=ntime;
    fprintf('Simulation stopped'); 
    break 
 end 
 
  dlmwrite('XYZ.txt',[X(1),X(2),X(3)], 'newline', 'pc', '-append')
  Y = dlmread('XYZ.txt'); set(h,'XData',Y(:,1),'YData',Y(:,2),'ZData',Y(:,3)) 
 end
