close all
clear all
clc

%% 2D Trajectory tracking simulation

% Simulating the trajectory
t = 0:0.01:29.99; T = length(t);
wt = 10./(t+5);
r = 10; 
xt = r*cos(wt.*t); yt = r*sin(wt.*t);
theta = [xt;yt]; % reference trajectory
% figure; plot(xt,yt); axis equal

% RHGD
mu = 1; beta = 5;
eta = 1e-2; % GD step size
gamma = 1e-2; % OGD step size

r = [];

W = 100; % Prediction window size

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



figure;
plot(theta(1,:),theta(2,:), 'LineWidth', 2); hold on
plot(x(1,:),x(2,:), 'LineWidth', 2);
xlabel('x (m)'); ylabel('y (m)'); axis equal
legend('target','agent')
p=strcat('Trajectory Tracking: Strongly Convex W=',int2str(W));
title(p)



% figure;
% for i = 1:T
%     plot(theta(1,1:i),theta(2,1:i),'LineWidth',3); drawnow
%     plot(x(1,1:i),x(2,1:i),'LineWidth',3); drawnow
%     axis equal
%     pause(0.001)
% end

%{
%% Initialize video
myVideo = VideoWriter('videoSC_50'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

   
figure;
title('Strongly Convex Cost')
subplot(2,2,2)
% plot(t, [theta(1,:);x(1,:)], 'LineWidth', 2);
hh1(1) = animatedline(t(1), theta(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b', 'MaximumNumPoints',1);
hh1(2) = animatedline(t(1), x(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0], 'MaximumNumPoints',1);
hh1(3) = animatedline('LineWidth',2,'Color','b');
hh1(4) = animatedline('LineWidth',2,'Color',[0 .5 0]);

xlim([0,10]); ylim([-100,100]);
xlabel('time (sec)'); ylabel('x (m)');
title('X location')


subplot(2,2,4)
% plot(t, [theta(2,:);x(2,:)], 'LineWidth', 2);

hh2(1) = animatedline(t(1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b', 'MaximumNumPoints',1);
hh2(2) = animatedline(t(1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0], 'MaximumNumPoints',1);
hh2(3) = animatedline('LineWidth',2,'Color','b');
hh2(4) = animatedline('LineWidth',2,'Color',[0 .5 0]);

xlim([0,10]); ylim([-100,100]);
xlabel('time (sec)'); ylabel('y (m)');
title('Y location')

subplot(2,2,[1,3])
% plot(theta(1,:), theta(2,:), 'LineWidth', 2); hold on
% plot(x(1,:), x(2,:), 'LineWidth', 2);

hh3(1) = animatedline(theta(1,1), theta(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b', 'MaximumNumPoints',1);
hh3(2) = animatedline(x(1,1), x(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0], 'MaximumNumPoints',1);
hh3(3) = animatedline('LineWidth',2,'Color','b');
hh3(4) = animatedline('LineWidth',2,'Color',[0 .5 0]);

legend('target','agent','target trajectory','agent trajectory')
xlabel('x (m)'); ylabel('y (m)'); axis equal
xlim([-100,100]);ylim([-200,200]);
title('Trajectory Tracking')

tic;     % start timing
for id = 1:T
   % Update XData and YData
%    set(hh1(1), 'XData', t(id)           , 'YData', theta(1, id));
%    set(hh1(2), 'XData', t(id)           , 'YData', x(1, id));
   addpoints(hh1(1),t(id),theta(1,id))
   addpoints(hh1(2),t(id),x(1,id))
   addpoints(hh1(3),t(id),theta(1,id))
   addpoints(hh1(4),t(id),x(1,id))
     
%    set(hh2(1), 'XData', t(id)           , 'YData', theta(2, id));
%    set(hh2(2), 'XData', t(id)           , 'YData', x(2, id));
    
   addpoints(hh2(1),t(id),theta(2,id))
   addpoints(hh2(2),t(id),x(2,id))
   addpoints(hh2(3),t(id),theta(2,id))
   addpoints(hh2(4),t(id),x(2,id))

%    set(hh3(1), 'XData', theta(1,id)     , 'YData', theta(2, id));
%    set(hh3(2), 'XData', x(1,id)         , 'YData', x(2, id));
   
   addpoints(hh3(1),theta(1,id),theta(2,id))
   addpoints(hh3(2),x(1,id),x(2,id))
   addpoints(hh3(3),theta(1,id),theta(2,id))
   addpoints(hh3(4),x(1,id),x(2,id))
   
   drawnow;
   
   frame = getframe(gcf); %get frame
   writeVideo(myVideo, frame);
   
end
fprintf('Animation (Smart update): %0.2f sec\n', toc);
close(myVideo)


% figure
% plot(0:W,r,'LineWidth',3)
% xlabel('Prediction Window size (W)'); ylabel('Reg(RHGD)')
% title('For Strongly convex OCO problems')

%}

function reg = regretSC(x, theta, mu, beta)
% calculates regret of strongly convex tracking problem
reg = 0.5*mu*trace((x-theta)'*(x-theta)) + 0.5*beta*trace(diff(x')'*diff(x'));
end




   
   
   
   
   
   
   
   
   
   
   
   
   
   
   