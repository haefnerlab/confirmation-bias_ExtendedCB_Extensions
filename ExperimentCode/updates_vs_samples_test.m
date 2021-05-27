wts = zeros(5,11);
err = zeros(5,11);
for i=20:20:100
    params = Model.newModelParams();
    params.updates = i;
    params.sensory_info = 0.55;
    [weights, errors, ~] = Model.plotPK(params);
    wts(i,:)=weights;
    err(i,:)=errors;
end
% close all;
figure(10000);hold on
for i=20:20:100
    txt = ['updates = ',num2str(i)];
    errorbar(1:params.frames, wts(i,1:end-1), err(i,1:end-1),'DisplayName',txt,'LineWidth',2);
%     errorbar(params.frames+1, wts(i,end), err(i,end), '-r');
    xlabel('time');
    ylabel('weight');
    %     legend('20' ,'40' ,'60' ,'80' ,'100')
end
hold off
legend show
figure(10001);hold on
for i=20:20:100
    txt = ['updates = ',num2str(i)];
    plot(1:params.frames, wts(i,1:end-1)/sum(wts(i,1:end-1)),'DisplayName',txt,'LineWidth',2);
%     errorbar(params.frames+1, wts(i,end), err(i,end), '-r');
    xlabel('time');
    ylabel('weight');
    %     legend('20' ,'40' ,'60' ,'80' ,'100')
end
hold off
legend show

wts_s = zeros(7,11);
err_s = zeros(7,11);
s = [1 5 15 25 50 100 200];
for i=1:length(s)
    params = Model.newModelParams();
    params.samples = i;
    params.sensory_info = 0.55;
    params.sensory_info = 0.9;
    [weights, errors, ~] = Model.plotPK(params);
    wts_s(i,:)=weights;
    err_s(i,:)=errors;
end
figure(20000);hold on
for i=1:length(s)
    txt = ['samples = ',num2str(s(i))];
    errorbar(1:params.frames, wts_s(i,1:end-1), err_s(i,1:end-1),'DisplayName',txt,'LineWidth',2);
%     errorbar(params.frames+1, wts(i,end), err(i,end), '-r');
    xlabel('time');
    ylabel('weight');
    %     legend('20' ,'40' ,'60' ,'80' ,'100')
end
hold off
legend show

figure(20001);hold on
for i=1:length(s)
    txt = ['samples = ',num2str(s(i))];
    plot(1:params.frames, wts_s(i,1:end-1)/sum(wts_s(i,1:end-1)), 'DisplayName',txt,'LineWidth',2);
%     errorbar(params.frames+1, wts(i,end), err(i,end), '-r');
    xlabel('time');
    ylabel('weight');
    %     legend('20' ,'40' ,'60' ,'80' ,'100')
end
hold off
legend show
