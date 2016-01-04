%% PUI and PIO for state and disturbance estimation
% Sam Nazari
% November 2015
% This illustrative example shows how to estimate the disturbance and
% state.
clear
clc
dy = []
A = [-1 1 1;1 0 0;0 1 -1];
B = [1;0;1];
C = [1 0 0;0 1 0]; % this is the orginal system that shafai proposed, but
%I am changing it to make the pair A1,C observable
%C = [1 0 1;0 1 0];
D = 0;
E = [1;0;0];

%% Construct the state space system
% Let us construct the system as a state space object for simulation.
sys = ss(A,B,C,D);
%% Step one: Check rank
% We check to see if the rank(CE) = rank(E) = 1
rank(C*E)
rank(E)
PCE = pinv(C*E);

%% Step two: Compute Observer Matricies
% Since we have no issues in the previous step, let us compute the observer
% matrices now:
H = E*inv((C*E)'*(C*E))*(C*E)'
T = eye(3)-H*C
A1 = T*A

%% Step three: check observability 
% We must check the system A1,C observability
rank(obsv(A1,C))
% Since rank(obsv(A1,C))<n, the pair (A1,C) is not observable.  Let us
% utilize the PBH test to see which eigenvalue of A1 is the culprit. 
l = eig(A1);
lam1=l(1),lam2=l(2),lam3=l(3)
plam1 = rank([lam1*eye(3)-A1;C])
plam2 = rank([lam2*eye(3)-A1;C])
plam3 = rank([lam3*eye(3)-A1;C])
% it can be seen that the eigenvalue that fails the PBH test is already in
% the left half plane (at -1).  Therefore, the pair (A1,C) is detectible
% and a UIO exists. 
%% Pole Placement 
% Pole placement is used to assign the observer poles
K1 = [3 0;-1 2;-1/2 0]
%% Finish observer design
% Finally, the F and K matricies are computed
F = A1-K1*C
K = K1 + F*H
eig(F)
%% Sim the system
%
% Set up simulation initial condiditons
x1_0 = 1;
x2_0 = 0;
x3_0 = 0;

d = 10;
TSIM = 50;
h=1/10;
sim('PUIOPPIOsim')

%% Estimate the disturbance
xHat = xHat';
yH = yH';
dy = zeros(2,length(yH));
CAxhat = CAxhat';
PCE = pinv(C*E);
dy(1,1:end-1)=diff(yH(1,:))
dy(2,1:end-1)=diff(yH(2,:))
dyH= dy./h;
dHat_ft = PCE*dyH;
dHat_st = PCE*CAxhat;

% yH2 = C*xHat';
% yHd2 = [diff(yH2(1,:))./h;diff(yH2(2,:))./h];
% CAxH= C*A*xHat';
% dHd = pinv(C*E)*(yHd2-CAxH(:,1:end-1));
% dHr = [diff(dist2wks)./h];
%% Plot the results
figure,
plot(tout,X(:,1),'b'),hold on
plot(tout,X(:,2),'g'),hold on
plot(tout,X(:,3),'r')
xlabel('Time in seconds'),ylabel('State'),ylim([-10,10]),grid on
title('Positive Dynamic System State'),grid on
legend('x_1','x_2','x_3')

figure,
plot(tout,xHat(1,:),'b'),hold on
plot(tout,xHat(2,:),'g'),hold on
plot(tout,xHat(3,:),'r')
xlabel('Time in seconds'),ylabel('State Estimate Error'),ylim([-10,10]),grid on
title('Positive Dynamic System State Estimates'),grid on
l=legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$');
set(l,'Interpreter','Latex');

figure
plot(tout,Xerr(:,1),'b'),hold on
plot(tout,Xerr(:,2),'g'),hold on
plot(tout,Xerr(:,3),'r')
xlabel('Time in seconds'),ylabel('State Estimate Error'),grid on
title('Positive Dynamic System and Positive UIO Observer State Estimation Error')
legend('x_1','x_2','x_3')

figure,
plot(tout,dist2wks,'b'),xlabel('Time in seconds'),grid on
title('Disturbance input to positive system')

figure,
plot(tout,dyH(1,:)),ylim([0,1000]),xlabel('Time in seconds'),grid on
tm1=title('$\dot{\hat{y}}(1)$ Estimate')
set(tm1,'Interpreter','Latex')

figure,
plot(tout,dyH(2,:)),ylim([0,1000]),xlabel('Time in seconds'),grid on
tm2=title('$\dot{\hat{y}}(2)$ Estimate')
set(tm2,'Interpreter','Latex')

figure,
plot(tout,-CAxhat(1,:)),xlabel('Time in seconds'),grid on,ylim([-10,10])
t1=title('The $CA\hat{x}(1)$ term');
set(t1,'Interpreter','Latex');

figure,
plot(tout,-CAxhat(2,:)),xlabel('Time in seconds'),grid on,ylim([-10,10])
t2=title('The $CA\hat{x}(2)$ term');
set(t2,'Interpreter','Latex');

figure,
plot(dyH(1,:)-CAxhat(1,:)),ylim([-100,100]),xlabel('Time in seconds'),grid on,
t3 = title('Disturbance Estimate $\hat{d}$')
set(t3, 'Interpreter','Latex')

figure,
plot(dyH(2,1:end-1)-CAxhat(2,1:end-1)),ylim([-10,10]),xlabel('Time in seconds'),grid on,
t4 = title('The term $\dot{\hat{y}}(2)-CA\hat{x}(2)$')
set(t4, 'Interpreter','Latex')



