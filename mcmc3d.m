tbegin = cputime;
prior_lower = 1e-4;
prior_upper = 5;
sd = 0.1;
% prior = makedist("Uniform","lower",10^prior_lower,"upper",10^prior_upper);
prior = makedist("Uniform","lower",prior_lower,"upper",prior_upper);
sigma_delta = 1;
mu_delta = [15,25];
prior_delta = makedist("Normal","mu",sum(mu_delta)/2,"sigma",sqrt(20)*sigma_delta);
delta = sum(mu_delta)/2;
beta_prior = 0.5;
theta = [0.5,delta, mu_delta, sigma_delta];
beta = beta_prior;
numchains = 1;
% tiledlayout(numchains+1,1)
% hold on
%     beta = 0.5;
iterations = 2e3;
% iterations = 5;
acc = 0;
rej = 0;
llf = 4;
sigsq = 0.3464;
% fprintf("d1 = %2.4f,d2 = %2.4f",distance(theta,0.05),distance(theta,0.005))
[path,dz] = simulate(theta);
db = distance(path, dz);
ll = zeros(iterations+1,1);
betas = zeros(iterations,1);
deltas = zeros(iterations,1);
ll(1) = llf/2 -(db^2/(2*sigsq));
% h = plot(beta,delta,'.');hold on;
for i = 1:iterations
    fprintf("%d ",i);
    betas(i) = beta;
    deltas(i) = delta;
    ll(i+1) = ll(i);
    while 1
    betanew = random("Normal",(beta), sd);
    if (betanew > prior_lower) && (betanew < prior_upper)
%         fprintf('%d Beta out of bounds Beta = %4.4f Delta = %2.4f Likelihood=%4.8f Alpha=%2.4f\n', i, beta, delta, ll(i+1), alpha);
        break
    end
    end
%     betanew = 10^betanew;
    while 1
    deltanew = random("Normal", delta, sigma_delta);
    if (deltanew > 5) && (deltanew <25)
        break
    end
    end
    theta(1:2) = [betanew,deltanew];
    [path, disp] = simulate(theta);
    ll(i+1) = loglikelihood(distance(path,disp),sigsq, llf);
    likelihood_factor = exp(ll(i+1)-ll(i));
    alpha = (pdf(prior,betanew)*pdf(prior_delta,deltanew)/(pdf(prior,beta)*pdf(prior_delta,delta)))*likelihood_factor;
    if acceptor(alpha) == 1
        
        beta = betanew;
        betas(i) = betanew;
        deltas(i) = deltanew;
        delta = deltanew;
        acc = acc + 1; 
%         paths(1:length(path),:,acc) = path; 
%         fprintf('%d Accepted\n Theta = (%4.4f,%2.4f) Likelihood=(%4.8f,%2.4f,%2.4f)\n', i, beta, delta, ll(i+1), alpha, likelihood_factor);
    else
        rej = rej+1;
%         fprintf('%d Rejected Likelihood=(%4.8f,%2.4f,%2.4f)\n',i, ll(i+1), alpha, likelihood_factor);
%         fprintf('Current\n Theta = (%4.4f,%2.4f) LogLikelihood=%4.8f\n', beta, delta, ll(i));
        ll(i+1) = ll(i);
    end
    fprintf('Time taken: %7.4f seconds\n',cputime-tbegin);
     fprintf('Current Theta = (%4.4f,%2.4f) LogLikelihood=%4.8f\n', beta, delta, ll(i+1));
%     set(h,'XData',betas(1:i),'YData',deltas(1:i)) ;
%     drawnow;

    
    
end
fprintf('Time taken to run %d simulations: %7.4f seconds\n',iterations,cputime-tbegin);

% betascurr_burnout = betas(iterations/4:iterations);
% cl = 0.7/numchains;
% X = 1:iterations;
% nexttile
%     prior = fitdist(betascurr_burnout(:),"Beta");
% cl = 0.7;
% plot(X,betas(X),Color=[cl 1-cl cl*0.5])   
% histogram(betas(iterations/4:iterations),10, "Normalization","pdf")
% nexttile
% beta_burnout = zeros(numchains*(3*iterations/4+1),1);
% beta_burnout(1:iterations-iterations/4+1,1) = betas(iterations/4:iterations,1);
% for i =2:numchains
%     beta_burnout(1+(i-1)*(iterations-iterations/4+1):i*(iterations-iterations/4+1),1) = betas(iterations/4:iterations,i);
% end
% pd = fitdist(beta_burnout(:),"Beta")
% X = 0:0.001:1;
% plot(X,pdf(pd,X),Color=[0 0.22 0.37])
% hold on
% % nexttile
% for i = 1:numchains
%     cl = 0.7*i/numchains;
%     curr_betadist = fitdist(betas(iterations/4:iterations,i),"Beta");
% %     plot(0:0.001:1 ,binopdf(s(i),trials,0:0.001:1)*101,"Color",[cl 1-cl cl*0.5],"LineStyle","--")
%     plot(0:0.001:1 ,pdf(curr_betadist,X),"Color",[cl 1-cl cl*0.5],"LineStyle","--")
% end
% pd.mean
function Zs = loglikelihood(disp,sigsq, llf)
Zs = llf/2-((disp)^2/(2*sigsq));
end
% function Zs = likelihood(beta,sigsq)
% Zs = exp(-(distance(beta)^2/(2*sigsq)));
% end
function dis = distance(trajectory,tardisp)
    Z = diff(trajectory);
    Z = Z.^2;
    Z = sum(Z,2);
    Z = sqrt(Z);
    dis = (sum(Z) - tardisp);
end
function [trajectory,delX] = simulate(theta)
beta = theta(1);
delta = theta(2);
control = 0;
if beta*delta < 0.5
    control = 1;
end
% mu_delta = theta(3:4);
% sigma_delta = theta(5);
% sqrb = sqrt(10*beta);

syms x y z v1 v2 v3 
ntime=1000;

% ssf = max(0.02*min(1,(sigma_delta/(delta-mu_delta(1)))^2),0.025*min(1,(sigma_delta/(delta-mu_delta(2)))^2));
% ssf = min(sqrb,1)*max(ssf, 0.02*min(1, (sigma_delta/(delta-(sum(mu_delta)/2)))^2));
fprintf("Simulating Theta = (%2.4f,%2.4f)\n",beta,delta);
tstart = cputime;
X = [10 40 5 0.1 0.1 0.1]; tau=[45 5 45];
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
% h = plot3(X(1),X(2),X(3),'-');hold on

% Draw target
%  ht = plot3(tau(1),tau(2),tau(3),' . ');set(ht,'MarkerSize',100);
% hold on 

% Draw Obsatcle
% [X1,Y1,Z1] =  sphere(20); 
% hobs = surf(o(1)+X1*rad,o(2)+Y1*rad,o(3)+Z1*rad);set(hobs,'MarkerSize',1, 'FaceColor','r');

% Button to early stop
% hstop = uicontrol('Style','pushbutton','String','Stop', 'Position',[30 0 80 10],'callback','earlystop = 1;'); 
% earlystop = 0;

% axis([0 50 0 50 0 50]) 
% grid on 

% delete('XYZ.txt')
Y = zeros(ntime,3);
Y(1,:) = X(1:3);
 for t=2:ntime
     drawnow 
var={x,y,z,v1,v2,v3};
varval={X(1),X(2),X(3),X(4),X(5),X(6)};
% varval={X(1),X(2),X(3),X(4),X(5),X(6)};
X4 = zeros(6,4);
X4(:,1) = X(:);
kx = zeros(6,4);
stepsize = 0.1;
if control
   w = (subs(W,var,varval));
    if w < 4
        stepsize = 0.05;
    else
        stepsize = 0.1;
    end
end
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

 
%  if earlystop
%     t=ntime;
%     fprintf('Simulation stopped'); 
%     break 
%  end
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
%   set(h,'XData',Y(1:t,1),'YData',Y(1:t,2),'ZData',Y(1:t,3)) 
 end
 for k = 1:3
 Y(t+1,k) = tau(k);
 end
%  delete(h);
trajectory = Y(1:t+1,:);
fprintf('Time taken for simulation: %2.4f\n',cputime-tstart);end

function z = acceptor(x)
    t = random("Uniform",0,1);
    if t < x
        z = 1;
    else
        z = 0;
    end
end
