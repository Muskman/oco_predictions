%% tt all comparison
close all
clear all
clc

%% 2D Trajectory tracking simulation
% Simulating the trajectory
t = 0:0.01:9.99; T = length(t);
wt = 10./(t+1);
r = 100; 
xt = r*cos(wt.*t); yt = r*sin(wt.*t);
theta = [xt;yt]; % reference trajectory
% figure; plot(xt,yt); axis equal

%% Strongly convex case

% RHGD
mu = 1; beta = 5;
eta = 1e-2; % GD step size
gamma = 1e-2; % OGD step size

r = [];

for W = 0:5:500
    x = zeros(size(theta)); x(:,1) = theta(:,1);
    for i = 1:T

       % initialization using OGD
       if i==1
           for j=2:W+1
               x(:,j) = x(:,j-1) - mu*gamma*(x(:,j-1)-theta(:,j-1));
           end
       else
           if i+W<=T
               x(:,i+W) = x(:,i+W-1) - mu*gamma*(x(:,i+W-1)-theta(:,i+W-1));
           end
       end

       % GD updates
       if i>1
           for j = min(i+W-1,T):-1:i
               if j~=T
                   g = mu*(x(:,j)-theta(:,j))+beta*(2*x(:,j)-x(:,j-1)-x(:,j+1));
               else
                   g = mu*(x(:,j)-theta(:,j))+beta*(x(:,j)-x(:,j-1));
               end
               x(:,j) = x(:,j) - eta*g;
           end
       end
    end
r = [r; regretSC(x,theta,mu,beta)];    
end

%{ 

figure;
plot(theta(1,:),theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:),x(2,:), 'LineWidth', 2);
xlabel('x (m)'); ylabel('y (m)'); axis equal
legend('target','agent')



% figure;
% for i = 1:T
%     plot(theta(1,1:i),theta(2,1:i),'LineWidth',3); drawnow
%     plot(x(1,1:i),x(2,1:i),'LineWidth',3); drawnow
%     axis equal
%     pause(0.001)
% end
   
   
figure;
subplot(2,2,2)
plot(t, [theta(1,:);x(1,:)], 'LineWidth', 2);
hh1(1) = line(t(1), theta(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(2) = line(t(1), x(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('x (m)');


subplot(2,2,4)
plot(t, [theta(2,:);x(2,:)], 'LineWidth', 2);
hh2(1) = line(t(1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh2(2) = line(t(1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('y (m)');

subplot(2,2,[1,3])
plot(theta(1,:), theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:), x(2,:), 'LineWidth', 2);
hh3(1) = line(theta(1,1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh3(2) = line(x(1,1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('x (m)'); ylabel('y (m)'); axis equal


tic;     % start timing
for id = 1:T
   % Update XData and YData
   set(hh1(1), 'XData', t(id)           , 'YData', theta(1, id));
   set(hh1(2), 'XData', t(id)           , 'YData', x(1, id));
   set(hh2(1), 'XData', t(id)           , 'YData', theta(2, id));
   set(hh2(2), 'XData', t(id)           , 'YData', x(2, id));
   set(hh3(1), 'XData', theta(1,id)     , 'YData', theta(2, id));
   set(hh3(2), 'XData', x(1,id)         , 'YData', x(2, id));

   drawnow;
end
fprintf('Animation (Smart update): %0.2f sec\n', toc);

%}


%% Quadratic Growth Case

% RHGD
mu = 1; beta = 5;
eta = 1e-2; % GD step size
gamma = 1e-2; % OGD step size

% huber function parameter
a = sqrt(0.5)*[5; 5];
rQG = [];

for W=0:5:500
    x = zeros(size(theta)); x(:,1) = theta(:,1);
    for i = 1:T
       % initialization using OGD
       if i==1
           for j=2:W+1
               if norm(x(:,j-1)-theta(:,j-1))>=norm(a)           
                   x(:,j) = x(:,j-1) - mu*gamma*(x(:,j-1)-theta(:,j-1)-a*sign(a'*(x(:,j-1)-theta(:,j-1))));
               else
                   x(:,j) = x(:,j-1);
               end
           end
       else
           if i+W<=T
               if norm(x(:,i+W-1)-theta(:,i+W-1))>=norm(a)
                   x(:,i+W) = x(:,i+W-1) - mu*gamma*(x(:,i+W-1)-theta(:,i+W-1)-a*sign(a'*(x(:,i+W-1)-theta(:,i+W-1))));
               else
                   x(:,i+W) = x(:,i+W-1);
               end
           end
       end

       % GD updates
       if i>1
           for j = min(i+W-1,T):-1:i
               if j~=T
                   if norm(x(:,j)-theta(:,j))>=norm(a)
                       g = mu*(x(:,j)-theta(:,j)-a*sign(a'*(x(:,j)-theta(:,j))))+beta*(2*x(:,j)-x(:,j-1)-x(:,j+1));
                   else
                       g = beta*(2*x(:,j)-x(:,j-1)-x(:,j+1));
                   end
               else
                   g = beta*(x(:,j)-x(:,j-1));
               end
               x(:,j) = x(:,j) - eta*g;
           end
       end
    end
    rQG = [rQG; regretQG(x,theta,beta,a)];
end


%{

figure;
plot(theta(1,:),theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:),x(2,:), 'LineWidth', 2);
xlabel('x (m)'); ylabel('y (m)'); axis equal
legend('target','agent')

% figure;
% for i = 1:T
%     plot(theta(1,1:i),theta(2,1:i),'LineWidth',3); drawnow
%     plot(x(1,1:i),x(2,1:i),'LineWidth',3); drawnow
%     axis equal
%     pause(0.001)
% end
   
   
figure;
subplot(2,2,2)
plot(t, [theta(1,:);x(1,:)], 'LineWidth', 2);
hh1(1) = line(t(1), theta(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(2) = line(t(1), x(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('x (m)');


subplot(2,2,4)
plot(t, [theta(2,:);x(2,:)], 'LineWidth', 2);
hh2(1) = line(t(1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh2(2) = line(t(1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('y (m)');

subplot(2,2,[1,3])
plot(theta(1,:), theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:), x(2,:), 'LineWidth', 2);
hh3(1) = line(theta(1,1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh3(2) = line(x(1,1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('x (m)'); ylabel('y (m)'); axis equal


tic;     % start timing
for id = 1:T
   % Update XData and YData
   set(hh1(1), 'XData', t(id)           , 'YData', theta(1, id));
   set(hh1(2), 'XData', t(id)           , 'YData', x(1, id));
   set(hh2(1), 'XData', t(id)           , 'YData', theta(2, id));
   set(hh2(2), 'XData', t(id)           , 'YData', x(2, id));
   set(hh3(1), 'XData', theta(1,id)     , 'YData', theta(2, id));
   set(hh3(2), 'XData', x(1,id)         , 'YData', x(2, id));

   drawnow;
end
fprintf('Animation (Smart update): %0.2f sec\n', toc);

%}


%% Convex case
% RHGD
mu = 1; beta = 5;
eta = 1e-2; % GD step size
gamma = 1e-2; % OGD step size

% huber function parameter
del = 5;
rConv = [];

for W=0:5:500
    x = zeros(size(theta)); x(:,1) = theta(:,1);
    for i = 1:T

       % initialization using OGD
       if i==1
           for j=2:W+1
               if norm(x(:,j-1)-theta(:,j-1))<=del           
                   x(:,j) = x(:,j-1) - mu*gamma*(x(:,j-1)-theta(:,j-1));
               else
                   x(:,j) = x(:,j-1) - sign(x(:,j-1)-theta(:,j-1))*del*gamma;
               end
           end
       else
           if i+W<=T
               if norm(x(:,i+W-1)-theta(:,i+W-1))<=del
                   x(:,i+W) = x(:,i+W-1) - mu*gamma*(x(:,i+W-1)-theta(:,i+W-1));
               else
                   x(:,i+W) = x(:,i+W-1) - sign(x(:,i+W-1)-theta(:,i+W-1))*del*gamma;
               end
           end
       end

       % GD updates
       if i>1
           for j = min(i+W-1,T):-1:i
               if j~=T
                   if norm(x(:,j)-theta(:,j))<=del
                       g = mu*(x(:,j)-theta(:,j))+beta*(2*x(:,j)-x(:,j-1)-x(:,j+1));
                   else
                       g = del*sign(x(:,j)-theta(:,j))+beta*(2*x(:,j)-x(:,j-1)-x(:,j+1));
                   end
               else
                   g = mu*(x(:,j)-theta(:,j))+beta*(x(:,j)-x(:,j-1));
               end
               x(:,j) = x(:,j) - eta*g;
           end
       end
    end
    rConv = [rConv; regretConv(x,theta,beta,del,mu)];
end
%{

figure;
plot(theta(1,:),theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:),x(2,:), 'LineWidth', 2);
xlabel('x (m)'); ylabel('y (m)'); axis equal
legend('target','agent')

% figure;
% for i = 1:T
%     plot(theta(1,1:i),theta(2,1:i),'LineWidth',3); drawnow
%     plot(x(1,1:i),x(2,1:i),'LineWidth',3); drawnow
%     axis equal
%     pause(0.001)
% end
   
   
figure;
subplot(2,2,2)
plot(t, [theta(1,:);x(1,:)], 'LineWidth', 2);
hh1(1) = line(t(1), theta(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(2) = line(t(1), x(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('x (m)');


subplot(2,2,4)
plot(t, [theta(2,:);x(2,:)], 'LineWidth', 2);
hh2(1) = line(t(1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh2(2) = line(t(1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('y (m)');

subplot(2,2,[1,3])
plot(theta(1,:), theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:), x(2,:), 'LineWidth', 2);
hh3(1) = line(theta(1,1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh3(2) = line(x(1,1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('x (m)'); ylabel('y (m)'); axis equal


tic;     % start timing
for id = 1:T
   % Update XData and YData
   set(hh1(1), 'XData', t(id)           , 'YData', theta(1, id));
   set(hh1(2), 'XData', t(id)           , 'YData', x(1, id));
   set(hh2(1), 'XData', t(id)           , 'YData', theta(2, id));
   set(hh2(2), 'XData', t(id)           , 'YData', x(2, id));
   set(hh3(1), 'XData', theta(1,id)     , 'YData', theta(2, id));
   set(hh3(2), 'XData', x(1,id)         , 'YData', x(2, id));

   drawnow;
end
fprintf('Animation (Smart update): %0.2f sec\n', toc);


%}


% Plot all together

figure
plot(0:5:W,r,'LineWidth',3); hold on
plot(0:5:W,rQG,'LineWidth',3)
plot(0:5:W,rConv,'LineWidth',3)
legend('Strongly Convex','Quad Growth','Convex')
xlabel('Prediction Window size (W)'); ylabel('Reg(RHGD)')
title('Regret variation with prediction window for OCO problems')



%% cost functions for various cases
function reg = regretSC(x, theta, mu, beta)
% calculates regret of strongly convex tracking problem
reg = 0.5*mu*trace((x-theta)'*(x-theta)) + 0.5*beta*trace(diff(x')'*diff(x'));
end

function reg = regretQG(x,theta,beta,a)
% regret for QG problems
reg = 0;
T = size(theta,2);
for i=1:T
    n = norm(x(:,i)-theta(:,i));
    if n>=norm(a)
        reg = reg + 0.5*(n^2+norm(a)^2-2*abs(a'*(x(:,i)-theta(:,i))));
    end
end
reg = reg + 0.5*beta*trace(diff(x')'*diff(x'));
end

function reg = regretConv(x,theta,beta,del,mu)
% regret for convex problems
reg = 0;
T = size(theta,2);
for i=1:T
    if norm(x(:,i)-theta(:,i))<del
        reg = reg+0.5*mu*norm(x(:,i)-theta(:,i))^2;
    else
        reg = reg+del*(norm(x(:,i)-theta(:,i))-del*0.5);
    end
end
reg = reg + 0.5*beta*trace(diff(x')'*diff(x'));
end



   
   
   
   
   
   
   
   
   
   
   
   
   
   
   