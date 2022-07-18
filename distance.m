% beta = theta(1);
% delta = theta(2);
beta = 0;
delta = 10;
% mu_delta = theta(3:4);
% sigma_delta = theta(5);
% sqrb = sqrt(10*beta);

syms x y z v1 v2 v3 
ntime=1000000;

% ssf = max(0.02*min(1,(sigma_delta/(delta-mu_delta(1)))^2),0.025*min(1,(sigma_delta/(delta-mu_delta(2)))^2));
% ssf = min(sqrb,1)*max(ssf, 0.02*min(1, (sigma_delta/(delta-(sum(mu_delta)/2)))^2));
fprintf("Theta = (%2.4f,%2.4f)\n",beta,delta);
tstart = cputime;
X = [10 40 5 0.1 0.1 0.1]; tau=[45 5 45];
% X = [10 40 5 0.1 0.1 0.1]; tau=[45 5 45];
delX = sqrt((X(1)-tau(1))^2  + (X(2)-tau(2))^2 + (X(3)-tau(3))^2 );
% delX = (X(1)-tau(1))^2  + (X(2)-tau(2))^2 + (X(3)-tau(3))^2;
o=[25 25 20]; rad=3;
delta1=delta;delta2=delta;delta3=delta; 

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
 ht = plot3(tau(1),tau(2),tau(3),' . ');set(ht,'MarkerSize',100);
hold on 

% Draw Obsatcle
[X1,Y1,Z1] =  sphere(20); 
 hobs = surf(o(1)+X1*rad,o(2)+Y1*rad,o(3)+Z1*rad);set(hobs,'MarkerSize',1, 'FaceColor','r');

% Button to early stop
% hstop = uicontrol('Style','pushbutton','String','Stop', 'Position',[30 0 80 10],'callback','earlystop = 1;'); 
earlystop = 0;

axis([0 50 0 50 0 50]) 
grid on 

% delete('XYZ.txt')
Y = zeros(ntime,3);
 for t=1:ntime
     drawnow 
var={x,y,z,v1,v2,v3};
% varval={X(1),X(2),X(3),X(4),X(5),X(6)};
% varval={X(1),X(2),X(3),X(4),X(5),X(6)};
X4 = zeros(6,4);
X4(:,1) = X(:);
kx = zeros(6,4);
for J = 1:4
    %     stepsize = min(1,12.5/delta)*min(1, delta/10)*0.05*min(1 , 0.02*exp(w/16)+exp(0.4*(w-16))); 
    stepsize = 0.05;
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
if (X(1)-tau(1))^2  + (X(2)-tau(2))^2 + (X(3)-tau(3))^2 < 400
%     fprintf(" Done");
    Y(t,1) = X(1);
    Y(t,2) = X(2);
    Y(t,3) = X(3);
    break
end
 Y(t,1) = X(1);
 Y(t,2) = X(2);
 Y(t,3) = X(3);
 
%   dlmwrite('XYZ.txt',[X(1),X(2),X(3)], 'newline', 'pc', '-append')
%   Y = dlmread('XYZ.txt');
  set(h,'XData',Y(1:t,1),'YData',Y(1:t,2),'ZData',Y(1:t,3)) 
 end
 for k = 1:3
 Y(t+1,k) = tau(k);
 end
%  delete(h);
Z = diff(Y(1:t+1,:));
Z = Z.^2;
Z = sum(Z,2);
Z = sqrt(Z);
Zs = (sum(Z) - delX);
fprintf('Time taken for simulation: %2.4f\n',cputime-tstart);
