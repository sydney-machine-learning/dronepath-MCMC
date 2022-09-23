tStart = cputime;
syms x y z v1 v2 v3 
ntime=1000000;stepsize = 0.1; 
X = [10 40 5 0.1 0.1 0.1]; tau=[45 5 45];

o=[25 25 25]; rad=4;
o2 = [15 35 10];
delta = 24.90;
delta1=delta;delta2=delta;delta3=delta; beta=0.1773; 

W=0.5*((x-o(1))^2+(y-o(2))^2+(z-o(3))^2-(rad)^2);
W2=0.5*((x-o2(1))^2+(y-o2(2))^2+(z-o2(3))^2-(rad)^2);
V=0.5*((x-tau(1))^2+(y-tau(2))^2+(z-tau(3))^2+v1^2+v2^2+v3^2); 
F=0.5*((x-tau(1))^2+(y-tau(2))^2+(z-tau(3))^2);

L=V+F*beta*(1/W+1/W2);
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
 [X2,Y2,Z2] =  sphere(20); 
 hobs = surf(o2(1)+X2*rad,o2(2)+Y2*rad,o2(3)+Z2*rad);set(hobs,'MarkerSize',1, 'FaceColor','b');

% Button to early stop
hstop = uicontrol('Style','pushbutton','String','Stop', 'Position',[30 0 80 10],'callback','earlystop = 1;'); 
earlystop = 0;

axis([0 50 0 50 0 50]) 
grid on 

delete('XYZ.txt');
% stepsize = 0.05;

 for t=1:ntime
     drawnow 
    var={x,y,z,v1,v2,v3};
    X4 = zeros(6,4);
    X4(:,1) = X(:);
    kx = zeros(6,4);
%     w = (subs(W,var,varval));
%     if w < 4
%         stepsize = 0.01;%*min(1,sqrt(0.04*beta*delta));
%     else
%         stepsize = 0.1;
%     end
        
for J = 1:4
    %     stepsize = min(1,12.5/delta)*min(1, delta/10)*0.05*min(1 , 0.02*exp(w/16)+exp(0.4*(w-16))); 
    
    % ODEs
    varval={X4(1,J),X4(2,J),X4(3,J),X4(4,J),X4(5,J),X4(6,J)};
    FX(1)=X(4);
    FX(2)=X(5);
    FX(3)=X(6);
    FX(4)=subs(sigma1,var,varval); 
    FX(5)=subs(sigma2,var,varval);
    FX(6)=subs(sigma3,var,varval);
    %     w = (subs(W,var,varval));
    % RK4 to numerically solve the ODE
    for I = 1:length(X)
        kx(I,J) = stepsize * FX(I);
        if J < 4
            X4(I,J+1) = X(I) + kx (I,J)*(0.5)*(1+floor(J/3));
        else
            X(I) = X(I) + (kx(I,1) + 2 * (kx(I,2) + kx(I,3)) + kx(I,4)) / 6;
        end
    end
end

 
 if earlystop
    t=ntime;
    fprintf('Simulation stopped'); 
    break 
 end 
 
  dlmwrite('XYZ.txt',[X(1),X(2),X(3)], 'newline', 'pc', '-append')
  Y = dlmread('XYZ.txt'); set(h,'XData',Y(:,1),'YData',Y(:,2),'ZData',Y(:,3)) 
 end
tEnd = cputime-tStart;