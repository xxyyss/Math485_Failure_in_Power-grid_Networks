%{
Test Functions for rn. 
    power_flow_f looks good. 
    
rn this is scrap to see if I can start figuring this out.
%}
clear all;
import const.*

a = 13;
b = 8;

idx = idxVal;
testCase = 'case14';
data = DataClass(testCase);


%Initilize Omega
omega_dot = zeros(1,data.nodes(3));
delta_dot = zeros(1,data.nodes(1));
neu_dot   = zeros(data.nodes(3),data.nodes(3));
neu   = zeros(data.nodes(3),data.nodes(3));
for i = 1:data.nodes(2)
     mi = data.network_data.branch(i,idx.FROM_BUS);
     mj = data.network_data.branch(i,idx.TO_BUS);
     neu(mi,mj) = .98;
end
neu(1,4) = 0;
neu(8,2) = 0;
%neu(3,6) = 0;
omega = zeros(1,data.nodes(3));

for i = 1:data.nodes(3)
    if ismember(i,data.gen)
        omega(i)  = 0;
    end
end
delta (:) = zeros(1,data.nodes(3));%randn(1,data.nodes(3))./10000;
delta(1)  = 0;

%
h = 0.001; %Step Size

%
fun = Functions();

for cnt = 1:10000
%k1 = fun.update_omega(data, omega, delta, h);
%k2 = h*fun.update_omega(data, omega + k1/2, delta, h);
%k3 = h*fun.update_omega(data, omega + k2/2, delta, h);
%k4 = h*fun.update_omega(data, omega + k3, delta, h);

    omega_dot = fun.update_omega(data, omega, delta, h);
    delta_dot = fun.update_delta(data, omega, delta, neu, h);
    neu_dot   = fun.update_neu( data, idx, a, b, delta, neu, h);
    omega = omega_dot*h + omega;
    delta = delta_dot*h + delta;
    neu   = neu_dot*h   + neu;
    neu
    omega
    delta
    cnt
    pause(0.1)
end 
