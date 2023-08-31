clear all
close all

files = {'08_29_2023_16_52_37','08_29_2023_16_52_06','08_29_2023_16_48_37','08_29_2023_16_48_20'};
titles = {'Coupled Stiffness'};
% 
d = dir;
isub = [d(:).isdir]; %# returns logical vector
files = {d(isub).name}';
files(ismember(files,{'.','..','old model'})) = [];

conversion = 0.00269*4.1/2.3;

L1 = 66.68/1000;
beta1 = 0;
beta2 = -16.43*3.14/180; 
y01 = 0;
y02 = 1.17/1000;
z01 = 33.1/2000;
z02 = -33.1/2000;
z0 = [z01;z02];
y0 = [y01; y02];
y1 = [L1; 66.81]/1000;
L2 = (y1(2) - y0(2))./cos(beta2);

% files = files(length(files)-38:length(files)-0);
N = length(files);
for i = 1:N
    load(files{i}+"/qs.mat");
    q_data = q_data.' - 3.14;
    load(files{i}+"/taus.mat");
    
    load(files{i}+"/timestamps.mat");
    
    load(files{i}+"/q_desired.mat");
    q_data(2,:) = -q_data(2,:);
    q1 = q_data(1,end);
    q2 = q_data(2,end);
    tau(i) = abs(tau_data(end)*conversion);
    Fk_al(i) = -(((L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(y01 - y02) + (L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(z01 - z02) + (2*(L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(L1*cos(beta1)*(cos(q1) - sin(beta1)^2*(cos(q1) - 1)) - L2*cos(beta2)*(cos(q2) - sin(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)*sin(beta1)^2*(cos(q1) - 1) + L2*cos(beta2)*sin(beta2)^2*(cos(q2) - 1)))/3 + (2*(L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(L1*sin(beta1)*(cos(q1) - cos(beta1)^2*(cos(q1) - 1)) - L2*sin(beta2)*(cos(q2) - cos(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)^2*sin(beta1)*(cos(q1) - 1) + L2*cos(beta2)^2*sin(beta2)*(cos(q2) - 1)))/3 - (2*(L1*cos(beta1)^2*cos(q1) - L1*cos(q1)*sin(beta1)^2)*(L1*sin(q1)*cos(beta1)^2 - L2*sin(q2)*cos(beta2)^2 - L1*sin(q1)*sin(beta1)^2 + L2*sin(q2)*sin(beta2)^2))/3))/2;

    k_al(i) = tau(i)/abs(Fk_al(i));
    
    
    Fk_r(i,:) = Fk_rej([q1, q2], beta1, beta2, 1, y0, y1, z0, [L1, L2]);
    k_rej(i) = tau(i)/abs(Fk_r(i,1));


    Fk_ar(i,:) = Fk_area([q1, q2], beta1, beta2, 1, y0, y1, z0, [L1, L2]);
    k_ar(i) = tau(i)/abs(Fk_ar(i,1));

    Fk_lin(i) = abs(q1 - q2);


    % tau{i} = tau_data(2:end);
    t{i} = time_data;

    q_end(:,i) = [q1;q2];
    q_{i} = q_data;
end

small = abs(q_end(1,:)-q_end(2,:)) < median(abs(q_end(1,:)-q_end(2,:)));
% small = abs(tau) < median(tau);
small = ~small;
small(small == 0) = 1;
k_al_avg = mean(k_al);
% std(k_al)

k_rej_avg = mean(k_rej);
% std(k_rej)

% mean(k_ar)
% std(k_ar)
disp('Results')
X_r = [abs(Fk_r(small,1)), ones(sum(small),1)];
y = abs(tau(small).');
K = (X_r.'*X_r) \ X_r.' * y;
r_r = y - X_r * K;
S = sum(r_r.^2);
R2 = sum((X_r * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Rejection model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)

X_al = [abs(Fk_al(small))', ones(sum(small),1)];
y = abs(tau(small).');
K = (X_al.'*X_al) \ X_al.' * y;
r_al = y - X_al * K;
S = sum(r_al.^2);
R2 = sum((X_al * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Distance model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)

X_lin = [abs(Fk_lin(small))', ones(sum(small),1)];
K = (X_lin.'*X_lin) \ X_lin.' * y;
r_lin = y - X_lin * K;
S = sum(r_lin.^2);
R2 = sum((X_lin * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Linear model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)


X_ar = [abs(Fk_ar(:,1)), ones(N,1)];
y = abs(tau.');
K = (X_ar.'*X_ar) \ X_ar.' * y;
r_ar = y - X_ar * K;
S = sum(r_ar.^2);
R2 = sum((X_ar * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Area model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)

% q1_test = linspace(-pi/2,pi/2);
% for i = 1:length(q1_test)
%     q1 = q1_test(i);
%     q2 = q_end
%     Fk_al(i) = -(((L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(y01 - y02) + (L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(z01 - z02) + (2*(L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(L1*cos(beta1)*(cos(q1) - sin(beta1)^2*(cos(q1) - 1)) - L2*cos(beta2)*(cos(q2) - sin(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)*sin(beta1)^2*(cos(q1) - 1) + L2*cos(beta2)*sin(beta2)^2*(cos(q2) - 1)))/3 + (2*(L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(L1*sin(beta1)*(cos(q1) - cos(beta1)^2*(cos(q1) - 1)) - L2*sin(beta2)*(cos(q2) - cos(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)^2*sin(beta1)*(cos(q1) - 1) + L2*cos(beta2)^2*sin(beta2)*(cos(q2) - 1)))/3 - (2*(L1*cos(beta1)^2*cos(q1) - L1*cos(q1)*sin(beta1)^2)*(L1*sin(q1)*cos(beta1)^2 - L2*sin(q2)*cos(beta2)^2 - L1*sin(q1)*sin(beta1)^2 + L2*sin(q2)*sin(beta2)^2))/3))/2;
% 
scatter(q_end(1,:),tau)
% [x,flag,relres,iter,resvec,lsvec] = lsqr(X_al,y)
%%
% find(t{1} > t{2}(end))
clear
close all
files = {'08_31_2023_15_40_41','08_31_2023_15_40_59'};
titles = {'Stiffness Compensation', 'No Stiffness Compensation'};
N = length(files);
conversion = 0.00269*4.1/2.3;
for i = 1:N
    load(files{i}+"/qs.mat");
    q_data = (q_data.' - 3.14);
    load(files{i}+"/taus.mat");
    
    load(files{i}+"/timestamps.mat");
    
    load(files{i}+"/q_desired.mat");
    q_data(2,:) = -q_data(2,:);
    q_desired = q_desired - pi;
    q1 = q_data(1,end);
    q2 = q_data(2,end);
    tau(i) = tau_data(end)*conversion;
    q_{i} = q_data;
    t{i} = time_data;
    figure
    
    colororder('default')
    plot(time_data, q_data(:,2:end), 'LineWidth', 1.3)
    hold on
    set(gca,'ColorOrderIndex',1)
    plot(time_data, repmat(q_desired,1,length(time_data)),'--', 'LineWidth', 1.3)
    title(titles{i})
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    ylabel('Joint Angle (rad)')
    % legend(module_legends)
    
    figure
    
    colororder('default')
    plot(time_data, tau_data(2:end), 'LineWidth', 1.3)
    hold on
    title(titles{i})
    set(gca,'FontSize',14)
    % xlim([0 20])
    xlabel('Time (s)')
    ylabel('Torque (Nm)')
    % legend(module_legends)
    set(gcf,'Position',[697,360,818,505])
end

