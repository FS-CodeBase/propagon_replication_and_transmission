function h = plot_one_data(propagon_data,sampling_times)
% % FUNCTION: plot_one_data
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
%
% DESCRIPTION 
% % Plot propagon amplification data by sampling times.
% 
% INPUTS
% % propagon_data: cell array of propagon data 
% % sampling_times: sampling times of propagon data
%
% OUTPUT
% % h: figure handle

for tidx = 1:length(sampling_times)
    T = sampling_times(tidx);
    for data = 1:numel(propagon_data{1,tidx})
        h = plot(T,propagon_data{1,tidx}(data),'-ko',...
               'MarkerFaceColor',[0.75,.75,.75],...
               'MarkerSize',8,...
               'Linewidth',1); hold on
    end
end
shg
hold off
xticks(0:8)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('Time (hrs)','Interpreter','latex','FontSize',19) 
ylabel('Number of Aggregates','Interpreter','latex','FontSize',18)
axis([0 8 0 700])
set(gcf, 'Position', [100, 100, 450, 800])
grid on
%print(str,'-dpng')