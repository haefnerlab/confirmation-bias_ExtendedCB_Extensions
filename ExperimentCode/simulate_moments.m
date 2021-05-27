close all;clear all;clc;
pi_c = linspace(0.45,0.55,5);
s = linspace(-0.1,0.1,5);
x = linspace(-1,1,1000);
mu_x = 0.25;
sigma_s = 0.5;
sigma_x = 0.5;
cnt=1;
mean_p_x_given_pi_c_and_s = zeros(length(pi_c),length(s));
var_p_x_given_pi_c_and_s = zeros(length(pi_c),length(s));
for i=1:length(pi_c)
    for j=1:length(s)
        for k=1:length(x)
            disp(cnt);
            cnt = cnt + 1;
            p_pos(k) = (normpdf(s(j),x(k),sigma_s) * normpdf(x(k),mu_x,sigma_x));
            p_neg(k) = (normpdf(s(j),x(k),sigma_s) * normpdf(x(k),-mu_x,sigma_x));            
        end
        p_pos = p_pos/sum(p_pos);
        p_neg = p_neg/sum(p_neg);
        for k=1:length(x)
            p_x_given_pi_c_and_s(i,j,k) = pi_c(i) * p_pos(k) + (1 - pi_c(i)) * p_neg(k);
            mean_p_x_given_pi_c_and_s(i,j) = mean_p_x_given_pi_c_and_s(i,j) + x(k)*p_x_given_pi_c_and_s(i,j,k);
            var_p_x_given_pi_c_and_s(i,j) = var_p_x_given_pi_c_and_s(i,j) + x(k)*x(k)*p_x_given_pi_c_and_s(i,j,k);
        end
% 
%         mean_p_x_given_pi_c_and_s(i,j) = mean(squeeze(p_x_given_pi_c_and_s(i,j,:)));
%         var_p_x_given_pi_c_and_s(i,j) = var(squeeze(p_x_given_pi_c_and_s(i,j,:)));
    end
end

%%
figure()
for i=1:length(pi_c)
    plot(mean_p_x_given_pi_c_and_s(i,:),var_p_x_given_pi_c_and_s(i,:),'o-','LineWidth',2);
    hold on;
end
xlabel('E_p[x] across s values')
ylabel('E_p[x^2] across s values');
legend (num2str(pi_c'))