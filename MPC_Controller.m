% quadcopter dynamic parameters

m=1.587; %mass
g=9.81; %gravity
Ixx=0.0213; %x-axis moment of inertia
Iyy=0.02217; %y-axis moment of inertia
Izz=0.0282; %z-axis moment of inertia
b=4.0687e-7; %thrust coefficient
k=8.4367e-9; %drag coefficient
L=0.243; %distance between the rotor and the center of mass
max_speed=4720; %max motor speed
min_speed=3093; %min take off speed

% MPC controller parameters

n_states=6; %number of states
n_inputs=3; %number of inputs
n_outputs=3; %number of outputs
nu=3; %control horizon
ny=6; %prediction horizon
n_sim=20; %number of simulation
x=zeros(6,1); %initial state values
y=zeros(3,1); %initial output values
u=[0;0;0]; %initial control values
t=0:0.2:n_sim; %step for sinusoidal reference

%reference to be tracked
% sin and cos reference
r=[sin(t)+cos(3*t)/2;sin(t)+cos(3*t)/2;sin(t)+cos(2*t)/2+sin(3*t)/3];
% sin path reference
% r=[sin(t);sin(t);sin(t)];
% helix path reference
% r=[2*cos(t);2*sin(t);t/(2*pi)];
Rs=[eye(n_outputs);eye(n_outputs);eye(n_outputs);...
    eye(n_outputs);eye(n_outputs);eye(n_outputs)]; %reference adjusting matrix
dist=(rand(3,101)*2-1)*0.5; %disturbance

%state space matrices
A=[0 1 0 0 0 0; ...
    0 0 0 0 0 0; ...
    0 0 0 1 0 0; ...
    0 0 0 0 0 0; ...
    0 0 0 0 0 1; ...
    0 0 0 0 0 0];

B=[0 0 0; 1/Ixx 0 0; 0 0 0;0 1/Iyy 0;0 0 0;0 0 1/Izz];
C=[0 1 0 0 0 0; ...
    0 0 0 1 0 0; ...
    0 0 0 0 0 1];

D=zeros(3,3);

%continuous to discrete time convertion
ts=0.2; %sampling time
[Ad,Bd,Cd,Dd]=c2dm(A,B,C,D,ts);

%augmented model
[Aa,Ba,Ca]=augment_mimo(Ad,Bd,Cd,n_states,n_inputs,n_outputs);

%controllability
CO=[Ba,Aa*Ba,Aa^2*Ba,Aa^3*Ba,Aa^4*Ba,Aa^5*Ba];
controllability=rank(CO);
if controllability==n_states
    disp("system is controllable")
end

%prediction matrices

%output predictions

P=[Ca*Aa;Ca*Aa^2;Ca*Aa^3;Ca*Aa^4;Ca*Aa^5;Ca*Aa^6];

H=[Ca*Ba,zeros(nu),zeros(nu);...
    Ca*Aa*Ba,Ca*Ba,zeros(nu);...
    Ca*Aa^2*Ba,Ca*Aa*Ba,Ca*Ba;...
    Ca*Aa^3*Ba,Ca*Aa^2*Ba,Ca*Aa*Ba;...
    Ca*Aa^4*Ba,Ca*Aa^3*Ba,Ca*Aa^2*Ba;...
    Ca*Aa^5*Ba,Ca*Aa^4*Ba,Ca*Aa^3*Ba];

% contraints
u1=L*b*(max_speed^2-min_speed^2);
u2=u1;
u3=k*(max_speed^2+max_speed^2-min_speed^2-min_speed^2);
u_max=[u1;u2;u3];
u_min=-u_max;
Dumax=0.6*u_max; %rate of input changes

%constraint matrices
[CC,dd,dupast]=constraints_mimo(Dumax,u_max,u_min,n_inputs,nu);

%cost functions

% weights on control inputs

W=[7.5e-2 0 0 0 0 0 0 0 0; ...
    0 7.5e-2 0 0 0 0 0 0 0; ...
    0 0 4.5e-2 0 0 0 0 0 0;...
    0 0 0 7.5e-2 0 0 0 0 0; ...
    0 0 0 0 7.5e-2 0 0 0 0; ...
    0 0 0 0 0 4.5e-2 0 0 0; ...
    0 0 0 0 0 0 7.5e-2 0 0; ...
    0 0 0 0 0 0 0 7.5e-2 0; ...
    0 0 0 0 0 0 0 0 4.5e-2];

E=2*(H'*H+W);

% simulation loop parameters
xf=zeros(size(Ba,1),1);
xh=zeros(size(Ba,1),1);
yh=y;

%simulation loop
for i = 1:100
    F=-2*H'*(Rs*r(:,i)-P*(xf));
    d= dd + dupast*u;
    [deltau,uncons]=Qphild(E,F,CC,d);
    
    %control horizon of 3 
    deltau=[deltau(1,1) deltau(2,1) deltau(3,1); ...
        deltau(4,1) deltau(5,1) deltau(6,1); ...
        deltau(7,1) deltau(8,1) deltau(9,1)];
    
    deltaU=deltau(1,:);
    u=u+deltaU'; %input
    
    %model
    xh(:,i+1)=Aa*xh(:,i)+Ba*deltaU';
    yh(:,i)=Ca*xh(:,i+1)+dist(:,i);
    xf=xh(:,i+1);    
end

% plots

 i=1:100;
    
% output1
subplot(3,1,1)
plot(i,r(1,i),'r',i,yh(1,i),'b')
legend('reference','output');
xlabel('time');
ylabel('phi')

% output2
subplot(3,1,2)
plot(i,r(2,i),'r',i,yh(2,i),'b')
legend('reference','output');
xlabel('time');
ylabel('theta')

% output3
subplot(3,1,3)
plot(i,r(3,i),'r',i,yh(3,i),'b')
legend('reference','output');
xlabel('time');
ylabel('psi')




