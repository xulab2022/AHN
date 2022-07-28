function fscope_beh_align_signal_mouse(fscopeMatPath,LogPath,FigPath,MatPath,video_fps,bins,bef_beh,beh_dur,filename)
% fscopeMatPath = '';
% filename = '';
% LogPath = '';
% FigPath = '';
% MatPath = '';
% video_fps = 25;
% bins = 20;  %ms
% bef_beh = 5; %s
% beh_dur = 10; %s

cd(fscopeMatPath)
info = load([filename,'.mat']);
fscope_signal = info.res.fscope_signal; %%%%
fscope_fps = info.res.fscope_fps;%%%%%

cd(LogPath)
info2 = load(['log',filename,'.mat']);%%%%
signal_start_frame = info2.markslog.framenum(2);
beh_start_frame = info2.markslog.framenum(3);
beh_end_frame = info2.markslog.framenum(4); 
action_code = info2.uc.bchar;
action_name = info2.uc.desc;

frame_idx = info2.rawlogframe;
action = info2.rawlog;
new_action = action;
action_idx = zeros(size(action,1),1);
for i = 1:size(action,1)
    if isspace(action{i})
        action_idx(i) = 1;
    end
end
act_idx = find(action_idx == 0);
new_action(action_idx == 1) = [];
if length(act_idx)>1
    tmp = diff(act_idx);
    for i = 1:size(tmp,1)
        if tmp(i) == 2
            act_idx2(i) = act_idx(i)+1;
        else act_idx2(i) = act_idx(i+1);
        end
    end
    frame(:,1) = frame_idx(act_idx);
    frame_tmp = frame_idx(act_idx2);
    frame_tmp(length(frame_tmp)+1) = frame_idx(end);
    frame(:,2) = frame_tmp;
    
elseif length(act_idx) == 1
    frame(:,1) = frame_idx(1);
    frame(:,2) = frame_idx(2);
else frame = [];
end


if ~isempty(find(diff(frame')<0, 1))
    error('The sequess of the action is not right')
end

indxx = [];
for i = 1:size(frame,1)
    if  frame(i,2) + beh_dur*video_fps > beh_end_frame
        indxx = [indxx,i];
    end
end
frame(indxx,:) = [];
new_action(indxx,:) = [];

new_frame = frame - signal_start_frame; 
fscope_frame = (new_frame)./video_fps.*fscope_fps+1; 
beh_start_fscope = (beh_start_frame - signal_start_frame)./video_fps.*fscope_fps + 1;
beh_end_fscope = (beh_end_frame - signal_start_frame)./video_fps.*fscope_fps + 1;
baseline = mean(fscope_signal(1 : 300001));%for 10min of+behavior

plot_start = round(fscope_frame(:,1) - bef_beh*fscope_fps); 
plot_end = round(fscope_frame(:,1) + beh_dur*fscope_fps); 
beh_plot = [plot_start,plot_end]; 

time_xaxis = -bef_beh:bins/1000:beh_dur;
for i = 1 : size(action_code,1);  %%%%%%%%%
    figure;
    subplot(2,1,1)
    indx = ismember(new_action,action_code{i});
    tmp = find(indx==1);
    if isempty(tmp) == 1;
        i = i+1;
    else if length(tmp) > 1;
        for j = 1:length(tmp)
            beh_plot_tmp(j,:) = (fscope_signal(beh_plot(tmp(j),1):round((fscope_fps*bins/1000)):beh_plot(tmp(j),2)) - baseline)./ baseline*100;    
        end
        mean_beh_plot = mean(beh_plot_tmp);
        std_beh_plot = std(beh_plot_tmp);
        sem_beh_plot = std_beh_plot./(sqrt(size(beh_plot_tmp,1)));
        L_sem = mean_beh_plot - sem_beh_plot;
        u_sem = mean_beh_plot + sem_beh_plot;
        for hhj = 1:length(time_xaxis)
            plot([time_xaxis(hhj) time_xaxis(hhj)],[L_sem(hhj) u_sem(hhj)],'color',[0.8 0.8 1]);
            hold on;
        end
        plot(time_xaxis,mean_beh_plot,'linewidth',1.5,'Color','b');
        hold on
        plot([0,0],ylim,'Color','r');
        set(gcf, 'position', get(0,'ScreenSize'));
        title([filename,'_',char(action_name{i}),],'interpreter','none','fontsize',15);
        cd(MatPath)
        res.mean_signal = mean_beh_plot;
        save([filename,'_',char(action_name{i}),'_'],'res');
    else if length(tmp) == 1;
            j = 1;
            beh_plot_tmp(1,:) = (fscope_signal(beh_plot(tmp(j),1):round((fscope_fps*bins/1000)):beh_plot(tmp(j),2)) - baseline)./ baseline*100;
            mean_beh_plot = beh_plot_tmp;
            plot(time_xaxis,mean_beh_plot,'linewidth',1.5,'Color','b');
            hold on
            plot([0,0],ylim,'Color','r');
            set(gcf, 'position', get(0,'ScreenSize'));
            title([filename,'_',char(action_name{i}),],'interpreter','none','fontsize',15);
            cd(MatPath)
            res.mean_signal = mean_beh_plot;
            save([filename,'_',char(action_name{i}),'_'],'res'); 
        end
        end
        subplot(2,1,2)
        imagesc(beh_plot_tmp)
        colormap('jet')
        caxis([-1 10]);
        colorbar('location','East')
        time_xaxis1 = 1:(1000/bins)*bef_beh:size(beh_plot_tmp,2);
        time_xaxis1_label = -bef_beh:5:beh_dur;
        set(gca,'XTick',time_xaxis1)
        set(gca,'XTickLabel',time_xaxis1_label);
        hold on
        plot([251,251],[0.5 (size(beh_plot_tmp,1)+0.5)],'r','linewidth',1.5);
        cd(FigPath)
        saveas(gcf,[filename,'_',char(action_name{i})],'bmp');   
   end
    indx = [];
    tmp = [];
    beh_plot_tmp = [];
    mean_beh_plot = [];
    std_beh_plot = [];
    sem_beh_plot = [];
    L_sem = [];
    u_sem = []; 
end

close all 
clear all




