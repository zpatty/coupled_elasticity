clear all
close all

files = {'08_29_2023_16_52_37','08_29_2023_16_52_06','08_29_2023_16_48_37','08_29_2023_16_48_20'};
titles = {'Coupled Stiffness'};
% 
d = dir;
isub = [d(:).isdir]; %# returns logical vector
files = {d(isub).name}';
files(ismember(files,{'.','..','old model','force control'})) = [];

conversion = 0.00269*4.1/2.3;
% conversion = 1;

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
l_vec = [L1;L2];
% 38
% files = files(length(files)-12-5:length(files)-4);
files = files(length(files)-38-5:length(files)-4);

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
    
    % Fk1 = calc_stiffness(q(1:2), beta1, beta2, k, y0, y1, z0, l_vec, flag);
    % Fk2 = calc_stiffness(flipud(q(1:2)), beta1, beta2, flipud(k), flipud(y0), flipud(y1), flipud(z0), flipud(l_vec), flag);
    % Fk = Fk1 - Fk2;
    Fk_r1(i,:) = Fk_rej(q1, q2, beta1, beta2, 1, y0, y1, z0, [L1, L2]);
    Fk_r2(i,:) = Fk_rej(q2, q1, beta2, beta1, 1, flipud(y0), flipud(y1), flipud(z0), [L2, L1]);
    Fk_r = Fk_r1 - Fk_r2;
    k_rej(i) = tau(i)/abs(Fk_r(i,1));


    Fk_ar(i,:) = Fk_area([q1, q2], beta1, beta2, 1, y0, y1, z0, [L1, L2]);
    k_ar(i) = tau(i)/abs(Fk_ar(i,1));

    Fk_lin(i) = abs(q1 - q2);
    
    

    % p1 = [L1*sin(beta1)^2*sin(q1) - L1*cos(beta1)^2*sin(q1); y01 + L1*cos(beta1)*(cos(q1) - sin(beta1)^2*(cos(q1) - 1)) - L1*cos(beta1)*sin(beta1)^2*(cos(q1) - 1); z01 + L1*sin(beta1)*(cos(q1) - cos(beta1)^2*(cos(q1) - 1)) - L1*cos(beta1)^2*sin(beta1)*(cos(q1) - 1)];
    % p2 = [L2*sin(beta2)^2*sin(q2) - L2*cos(beta2)^2*sin(q2); y02 + L2*cos(beta2)*(cos(q2) - sin(beta2)^2*(cos(q2) - 1)) - L2*cos(beta2)*sin(beta2)^2*(cos(q2) - 1); z02 + L2*sin(beta2)*(cos(q2) - cos(beta2)^2*(cos(q2) - 1)) - L2*cos(beta2)^2*sin(beta2)*(cos(q2) - 1)];
    % w = norm([y01 + L1*cos(beta1); z02 + L2*sin(beta2)] - [y02 + L2*cos(beta2); z01 + L1*sin(beta1)]);
    dtheta = abs(q1 - q2);
    if dtheta == 0
        lam = 0.0000000001;
    else
        w = 0.08/2;
        Px = L1*sin(dtheta);
        % Px = norm(p1 - p2)/2;
        alpha = atan2(Px,w);
        ell_1 = w/cos(alpha);
        lam = ell_1/w;
    end
    Fk_nh(i) = (lam - 1/lam^3);


    % tau{i} = tau_data(2:end);
    t{i} = time_data;

    q_end(:,i) = [q1;q2];
    q_{i} = q_data;
end

small = abs(q_end(1,:)-q_end(2,:)) < 1/3*median(abs(q_end(1,:)-q_end(2,:)));
% small = abs(q_end(1,:)-q_end(2,:))*180/pi < 45;

% small = abs(tau) < median(tau);
% small = ~small;
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
k_rej = K;

X_al = [abs(Fk_al(small))', ones(sum(small),1)];
y = abs(tau(small).');
K = (X_al.'*X_al) \ X_al.' * y;
r_al = y - X_al * K;
S = sum(r_al.^2);
R2 = sum((X_al * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Distance model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)
k_al = K;

X_lin = [abs(Fk_lin(small))', ones(sum(small),1)];
K = (X_lin.'*X_lin) \ X_lin.' * y;
r_lin = y - X_lin * K;
S = sum(r_lin.^2);
R2 = sum((X_lin * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['Linear model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)
k_lin = K;

% X_ar = [abs(Fk_ar(:,1)), ones(N,1)];
% y = abs(tau.');
% K = (X_ar.'*X_ar) \ X_ar.' * y;
% r_ar = y - X_ar * K;
% S = sum(r_ar.^2);
% R2 = sum((X_ar * K - mean(y)).^2)/sum((y - mean(y)).^2);
% s = ['Area model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
% disp(s)
% k_nh = K;

X_nh = [abs(Fk_nh(small))', ones(sum(small),1)];
K = (X_nh.'*X_nh) \ X_nh.' * y;
r_nh = y - X_nh * K;
S = sum(r_nh.^2);
R2 = sum((X_nh * K - mean(y)).^2)/sum((y - mean(y)).^2);
s = ['NeoH model  k = ', num2str(K'),', S = .', num2str(S), ' R^2 = ', num2str(R2)];
disp(s)
k_nh = K;

% q1_test = linspace(-pi/2,pi/2);
% for i = 1:length(q1_test)
%     q1 = q1_test(i);
%     q2 = q_end
%     Fk_al(i) = -(((L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(y01 - y02) + (L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(z01 - z02) + (2*(L1*cos(beta1)*(sin(q1) - sin(beta1)^2*sin(q1)) - L1*cos(beta1)*sin(beta1)^2*sin(q1))*(L1*cos(beta1)*(cos(q1) - sin(beta1)^2*(cos(q1) - 1)) - L2*cos(beta2)*(cos(q2) - sin(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)*sin(beta1)^2*(cos(q1) - 1) + L2*cos(beta2)*sin(beta2)^2*(cos(q2) - 1)))/3 + (2*(L1*sin(beta1)*(sin(q1) - cos(beta1)^2*sin(q1)) - L1*cos(beta1)^2*sin(beta1)*sin(q1))*(L1*sin(beta1)*(cos(q1) - cos(beta1)^2*(cos(q1) - 1)) - L2*sin(beta2)*(cos(q2) - cos(beta2)^2*(cos(q2) - 1)) - L1*cos(beta1)^2*sin(beta1)*(cos(q1) - 1) + L2*cos(beta2)^2*sin(beta2)*(cos(q2) - 1)))/3 - (2*(L1*cos(beta1)^2*cos(q1) - L1*cos(q1)*sin(beta1)^2)*(L1*sin(q1)*cos(beta1)^2 - L2*sin(q2)*cos(beta2)^2 - L1*sin(q1)*sin(beta1)^2 + L2*sin(q2)*sin(beta2)^2))/3))/2;
% 
scatter(q_end(2,end-13:end)-q_end(1,end-13:end),tau(end-13:end))
hold on
q = linspace(-0.26,0.6);
q2_plt = 0.76;
plot(q2_plt-q,abs(k_lin(1)*(q - q2_plt)) + k_lin(2),'LineWidth',1.3)
w = 0.08/2;
dtheta = abs(q - q2_plt);
Px = L1*sin(dtheta);
% Px = norm(p1 - p2)/2;
alpha = atan2(Px,w);
ell_1 = w./cos(alpha);
lam = ell_1/w;
plot(q2_plt-q,abs(k_nh(1)*(lam - 1./lam.^3))+k_nh(2),'LineWidth',1.3);
fkal = (((L1.*cos(beta1).*(sin(q) - sin(beta1)^2.*sin(q)) - L1.*cos(beta1).*sin(beta1)^2.*sin(q)).*(y01 - y02) + (L1.*sin(beta1).*(sin(q) - cos(beta1)^2.*sin(q)) - L1.*cos(beta1)^2.*sin(beta1).*sin(q)).*(z01 - z02) + (2.*(L1.*cos(beta1).*(sin(q) - sin(beta1)^2.*sin(q)) - L1.*cos(beta1).*sin(beta1)^2.*sin(q)).*(L1.*cos(beta1).*(cos(q) - sin(beta1)^2.*(cos(q) - 1)) - L2.*cos(beta2).*(cos(q2_plt) - sin(beta2)^2.*(cos(q2_plt) - 1)) - L1.*cos(beta1).*sin(beta1)^2.*(cos(q) - 1) + L2.*cos(beta2).*sin(beta2)^2.*(cos(q2_plt) - 1)))/3 + (2.*(L1.*sin(beta1).*(sin(q) - cos(beta1)^2.*sin(q)) - L1.*cos(beta1)^2.*sin(beta1).*sin(q)).*(L1.*sin(beta1).*(cos(q) - cos(beta1)^2.*(cos(q) - 1)) - L2.*sin(beta2).*(cos(q2_plt) - cos(beta2)^2.*(cos(q2_plt) - 1)) - L1.*cos(beta1)^2.*sin(beta1).*(cos(q) - 1) + L2.*cos(beta2)^2.*sin(beta2).*(cos(q2_plt) - 1)))/3 - (2.*(L1.*cos(beta1)^2.*cos(q) - L1.*cos(q).*sin(beta1)^2).*(L1.*sin(q).*cos(beta1)^2 - L2.*sin(q2_plt).*cos(beta2)^2 - L1.*sin(q).*sin(beta1)^2 + L2.*sin(q2_plt).*sin(beta2)^2))/3))/2;
plot(q2_plt-q,abs(k_al(1)*fkal) ,'LineWidth',1.3);
for i = 1:length(q)
    Fk_r1 = Fk_rej(q(i), q2_plt, beta1, beta2, 1, y0, y1, z0, [L1, L2]);
    Fk_r2 = Fk_rej(q2_plt, q(i), beta2, beta1, 1, flipud(y0), flipud(y1), flipud(z0), [L2, L1]);
    fkrej(i) = Fk_r1(1) - Fk_r2(1);
end
plot(q2_plt-q,abs(k_rej(1)*fkrej)+k_rej(2),'LineWidth',1.3);
legend({'', 'Linear', 'Neo-Hookean Shear', 'Distance', 'Rejection'})
xlabel('Angle between joints (rads)')
ylabel('Coupling Torque (Nm)')
set(gca,'FontSize',14)
% set(gca(), 'LooseInset', get(gca(), 'TightInset'));

% scatter(q_end(1,:),-Fk_r(:,1)*k_rej(1) + k_rej(2))
% [x,flag,relres,iter,resvec,lsvec] = lsqr(X_al,y)
%%
% find(t{1} > t{2}(end))
clear
close all
files = {'08_31_2023_15_40_41_comp','08_31_2023_15_40_59_no_comp'};
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
    plot(time_data, q_data(:,2:end), 'LineWidth', 2.3)
    hold on
    set(gca,'ColorOrderIndex',1)
    plot(time_data, repmat(q_desired,1,length(time_data)),'--', 'LineWidth', 2.2)
    % title(titles{i})
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    ylabel('Joint Angle (rad)')
    % set(gcf,'Position',[697,360,818,505])
    legend({'\theta_1', '\theta_2','\theta_1 ref'})
    % legend(module_legends)
    % xlim([0 6.6])
    figure
    
    colororder('default')
    plot(time_data, tau_data(2:end), 'LineWidth', 2.2)
    hold on
    % title(titles{i})
    set(gca,'FontSize',14)
    % xlim([0 20])
    xlabel('Time (s)')
    ylabel('Torque (Nm)')
    % xlim([0 6.6])
    % legend(module_legends)
    % set(gcf,'Position',[697,360,818,505])
end

%% Force Control

clear
close all
files = {'force control/09_05_2023_15_41_13_10', 'force control/09_05_2023_15_41_35_20', 'force control/09_05_2023_15_41_53_30', 'force control/09_05_2023_15_42_15_40'};
% titles = {'10', 'No Stiffness Compensation'};
N = length(files);
conv = 7.1073;
conversion = 0.00269*4.1/2.3/0.07;
% conversion = 1;
tau_d = [10; 20; 30; 40]*conversion;
figure
hold on
for i = 1:N
    % q_data = (q_data.' - 3.14);
    load(files{i}+"/taus.mat");
    load(files{i}+"/Fk.mat");
    load(files{i}+"/Fz.mat");
    
    load(files{i}+"/timestamps.mat");
    
    % q_data(2,:) = -q_data(2,:);
    % q_desired = q_desired - pi;
    % q1 = q_data(1,end);
    % q2 = q_data(2,end);
    % tau(i) = tau_data(end)*conversion;
    % q_{i} = q_data;
    
    t{i} = time_data;
    F_model{:,i} = F_model_data;
    F_ft{:,i} = F_fb_data;
    avg_force(i) = mean(F_fb_data);
    avg_torque(i) = mean(F_model_data);

    set(gca,'ColorOrderIndex',1)
    % figure
    plot(time_data(time_data < 6), conversion*F_model_data(time_data < 6), 'LineWidth', 2)
    % hold on
    plot(time_data(time_data < 6), conversion*smoothdata(F_fb_data(time_data < 6))*conv, 'LineWidth', 2)
    % colororder('default')
    set(gca,'ColorOrderIndex',2)
    plot(time_data(time_data < 6), repmat(tau_d(i),1,length(time_data(time_data < 6))), '--','LineWidth', 2)
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    ylabel('Force (N)')
    legend({'Estimated Force', 'Actual Force', 'Desired Force'})
    ylim([0 5])
        
end

figure
plot(avg_force,avg_torque);
conv = mean(avg_torque./avg_force);

%%
clear
close all
files = {'force control/09_05_2023_15_41_13_10', 'force control/09_05_2023_15_41_35_20', 'force control/09_05_2023_15_41_53_30', 'force control/09_05_2023_15_42_15_40'};
% titles = {'10', 'No Stiffness Compensation'};
N = length(files);
conv = 7.1073;
conversion = 0.00269*4.1/2.3/0.07;
% conversion = 1;
tau_d = [10; 20; 30; 40]*conversion;
figure
hold on
for i = 1:N
    % q_data = (q_data.' - 3.14);
    load(files{i}+"/taus.mat");
    load(files{i}+"/Fk.mat");
    load(files{i}+"/Fz.mat");
    
    load(files{i}+"/timestamps.mat");
    
    % q_data(2,:) = -q_data(2,:);
    % q_desired = q_desired - pi;
    % q1 = q_data(1,end);
    % q2 = q_data(2,end);
    % tau(i) = tau_data(end)*conversion;
    % q_{i} = q_data;
    t{i} = time_data;
    F_model{:,i} = F_model_data;
    F_ft{:,i} = F_fb_data;
    avg_force(i) = mean(F_fb_data);
    avg_torque(i) = mean(F_model_data);

    % set(gca,'ColorOrderIndex',1)
    % figure
    % plot(time_data, conversion*F_model_data(2:end), 'LineWidth', 2)
    % hold on
    plot(time_data(time_data < 6), conversion*smoothdata(F_fb_data(time_data < 6))*conv - tau_d(i), 'LineWidth', 2)
    % colororder('default')
    set(gca,'FontSize',14)
    xlabel('Time (s)')
    ylabel('Error (N)')
    ylim([-1 1.5])
    legend_str{i} = num2str(round(tau_d(i),1,'decimals')) + " N";
        
end
legend(legend_str)