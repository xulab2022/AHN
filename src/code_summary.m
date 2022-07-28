%% part1 for calculating each zone signal of open field and novel object
FscopeMatpath = '';  %fscopemat
MousepointMatpath = ''; % mousepointmat
LogPath = ''; % log
ExcelPath = ''; % excel
xlsname = 'test.xlsx';
video_fps = 25; 
zone = [20,30,40]; 
beh_time = 600;
filelist = sortdir(LogPath,'*.mat');
for i = 1:length(filelist);
    filename = filelist{i}(4:end-4);
    fscope_of_no_combine_during_record_zone_signal(FscopeMatpath,MousepointMatpath,LogPath,ExcelPath,xlsname,filename,video_fps,zone,beh_time)
end
%% part2 for aligning behavior with signal of each trial of each mouse
fscopeMatPath = ''; %fscopemat
LogPath = ''; %log
FigPath = ''; % mouse_fig
MatPath = ''; %mouse_mat
video_fps = 25;
bins = 20;  
bef_beh = 5; 
beh_dur = 10; 
filelist = sortdir(LogPath,'*.mat');
for i = 1:length(filelist);
    filename = filelist{i}(4:end-4);
    fscope_beh_align_signal_mouse(fscopeMatPath,LogPath,FigPath,MatPath,video_fps,bins,bef_beh,beh_dur,filename)
end
%% part3 for avearge the signal aligned with behavior of all trials from all mice
MatPath1 = ''; %vmouse_mat
MatPath2 = ''; % mice_mat
FigPath = '';  %mice_fig
bef_beh = 5;
beh_dur = 10;
bins = 20; 

cd(MatPath1)
file_list = sortdir(MatPath1,'*.mat');
for i = 1:length(file_list)
    indx = strfind(file_list{i},'_');
    virus{i,1} = file_list{i}(indx(2)+1 : (indx(3)-1)); 
    stimuli{i,1} = file_list{i}(indx(4)+1 : (indx(5)-1));
    end_beh{i,1} = file_list{i}((indx(end-1)):(end-4));
    indx = [];
end

viru = unique(virus);
stimu = unique(stimuli);
beh_type = unique(end_beh);

time_xaxis = -bef_beh:bins/1000:beh_dur; 

for i = 1:length(beh_type)   
    figure;
    filelist = sortdir(MatPath1,['*',beh_type{i},'.mat']); 
    if length(filelist) ~= 1;
        for j = 1:length(filelist)
            filename = (filelist{j}(1:end-4));
            cd(MatPath1)
            info = load([filename,'.mat']);
            beh_signal(j,:) = info.res.mean_signal;
        end
        ave_beh_signal = mean(beh_signal);
        std_beh_signal = std(beh_signal);
        sem_beh_signal = std_beh_signal./(sqrt(size(beh_signal,1)));
        L_sem = ave_beh_signal - sem_beh_signal;
        u_sem = ave_beh_signal + sem_beh_signal;
        for hhj = 1:length(time_xaxis)
            plot([time_xaxis(hhj) time_xaxis(hhj)],[L_sem(hhj) u_sem(hhj)],'color',[0.8 0.8 1]);
            hold on;
        end
    plot(time_xaxis,ave_beh_signal,'linewidth',1.5,'Color','b');
    else j = 1;
        filename = filelist{j}(1:end-4);
        cd(MatPath1)
        info = load([filename,'.mat']);
        beh_signal = info.res.mean_signal;
        ave_beh_signal = beh_signal;
        std_beh_signal = zeros(size(ave_beh_signal));
        sem_beh_signal = zeros(size(ave_beh_signal));
        L_sem = ave_beh_signal - sem_beh_signal;
        u_sem = ave_beh_signal + sem_beh_signal;
        plot(time_xaxis,ave_beh_signal,'linewidth',1.5,'Color','b');
    end
        hold on
        ylim([-2 10]);
        set(gca,'ytick',-2:4:10);
        plot([0,0],ylim,'Color','r');
        set(gcf, 'position', get(0,'ScreenSize'));
        title([cell2str(viru),'_',cell2str(stimu),beh_type{i},'average'],'interpreter','none','fontsize',15); 
        cd(MatPath2)
        res.beh_signal = beh_signal;
        res.ave_signal = ave_beh_signal;
        res.std_signal = std_beh_signal;
        res.sem_signal = sem_beh_signal;
        res.L_sem = L_sem;
        res.u_sem = u_sem;
        save([cell2str(viru),'_',cell2str(stimu),beh_type{i},'average'],'res');
        cd(FigPath)
        saveas(gcf,[cell2str(viru),'_',cell2str(stimu),beh_type{i},'average'],'bmp');
        print(gcf, '-dpdf', '-r600', [cell2str(viru),'_',cell2str(stimu),beh_type{i},'_average.pdf']);
        filelist = [];
        filename = [];
        info = [];
        beh_signal = [];  
        ave_beh_signal = [];
        std_beh_signal = [];
        sem_beh_signal = [];
        L_sem = [];
        u_sem = [];
end
close all
clear all
%% part4 for calculating each zone signal of epm
MatPath = ''; % fscopemat
MatPath2 = ''; % mousepoint_squart_mat
MatPath3 = ''; % epm_beh_info_mat
LogPath = ''; % log
ExcelPath = ''; % excel
xlsname = 'test.xlsx';
framerate = 25;

filelist = sortdir(LogPath,'*.mat');
for i = 1:length(filelist);
    filename = filelist{i}(4:end-4);
    fscope_epm_zone_signal(MatPath,MatPath2,MatPath3,LogPath,ExcelPath,xlsname,framerate,filename)
end
%% part5 Heatmap of epm signal
FscopeMatPath = '';  %fscopemat
MatPath = ''; %mat
LogPath = ''; %log
FigPath = ''; %figure
beh_time = 10; 
nbins = 5;
file_list = sortdir(LogPath,'*.mat');
color_caxis = [-5 8];
for i = 1:length(file_list)
    FileName = file_list{i}(4:end-4);
    fscope_align_sigal_EPM_heatmap(FscopeMatPath,MatPath,LogPath,FigPath,FileName,nbins,beh_time,color_caxis)
end
  

%% part6 beh_align_spike
[f1,p1] = uigetfile('*.mat'); 
[f2,p2] = uigetfile('*-wv.mat','MultiSelect','on'); 
FigPath = ''; % figure
fps = 30;
bins = 0.25; 
bef_beh = 5; 
post_beh = 5; 

info = load([p1,f1]);
beh_start = info.markslog.framenum(2);

beh_start_time = beh_start/fps;

 
action_code = info.uc.bchar;
action_name = info.uc.desc;

frame_idx = info.rawlogframe;
action = info.rawlog;
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

a_idx = strcmp(new_action,'a');
frame_appro_start = frame(a_idx,1)/fps;

for i = 1:length(frame_appro_start)
    range2(i,:) = (frame_appro_start(i,1) - bef_beh) : bins : (frame_appro_start(i,1) + post_beh);
end

for i = 1:length(f2)
    spike_mat = load([p2,cell2mat(f2(i))]);
    iTime = spike_mat.iTime;
    event_spike_bin = [];
    for h = 1:size(range2,1)
        for k = 1:size(range2,2)-1
            event_spike_bin_tmp = iTime(iTime > range2(h,k) & iTime <= range2(h,k+1));
            event_spike_bin = [event_spike_bin;event_spike_bin_tmp];
            beh_spike(h,k) = length(iTime(iTime > range2(h,k) & iTime <= range2(h,k+1)));
            beh_fr(h,k) = beh_spike(h,k)/bins;
        end
        event{h,:} = event_spike_bin;
    end
    mean_beh_fr = mean(beh_fr);
    figure
    subplot(2,1,1)
    for g = 1: size(event,1)
        x_indx_tmp = event{g,1} - frame_appro_start(g,1);
        for gg = 1:length(x_indx_tmp)
        y_indx1 = length(frame_appro_start)-g;
        y_indx2 = length(frame_appro_start)-g+1;
        plot([x_indx_tmp(gg) x_indx_tmp(gg)],[y_indx1 y_indx2],'color','k')
        hold on
        end
        plot([0 0],[0,length(frame_appro_start)],'color','r','linewidth',0.8)
        xlim([-bef_beh,post_beh])
        title(f2{i}(1:end-4),'interpreter','none','fontsize',15)
    end
    
    subplot(2,1,2)
    x_indx = -bef_beh+bins/2 : bins : post_beh-bins/2 ;
    bar(x_indx,mean_beh_fr,1,'grouped','Facecolor',[0 205/255 205/255])
    xlim([-bef_beh,post_beh])
    hold on
    plot([0 0],ylim,'color','r','linewidth',0.8);
    set(gcf, 'position', get(0,'ScreenSize'));
    cd(FigPath)
    saveas(gcf,f2{i}(1:end-4),'bmp');
   print(gcf,'-dpdf','-r600',[f2{i}(1:end-4),'.pdf']);
end
clear all
close all
disp('finish')
%% part7 for judging the neuron is excited or inhibited 
[f1,p1] = uigetfile('*.mat');
[f2,p2] = uigetfile('*-wv.mat','MultiSelect','on');
ExcelPath = ''; %excel
MatPath = ''; %mat
fps = 30;
bins = 0.25; 
beh_dur = 5;
bef_beh = 5; 
bef_beh2 = 0; 
beh_judge = 2;  
info = load([p1,f1]);
action_code = info.uc.bchar;
action_name = info.uc.desc;

frame_idx = info.rawlogframe;
action = info.rawlog;
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

a_idx = strcmp(new_action,'a');
frame_appro_start = frame(a_idx,1)/fps;

for i = 1:length(frame_appro_start)
    range(i,:) = (frame_appro_start(i,1) - bef_beh) : bins : (frame_appro_start(i,1) - bef_beh2);
end

for i = 1:length(frame_appro_start)
    range2(i,:) = frame_appro_start(i,1) : bins : (frame_appro_start(i,1) + beh_dur);
end

for i = 1:length(f2)
    spike_mat = load([p2,cell2mat(f2(i))]);
    iTime = spike_mat.iTime;
    %  phase1_itime = iTime(iTime <= (beh_time*60+1/fps));

    for j = 1:size(range,1)
        for t = 1:size(range,2)-1
            base_spike(j,t) = length(iTime(iTime > range(j,t) & iTime <= range(j,t+1)));
            base_fr(j,t) = base_spike(j,t)/bins;
        end
    end
    base_fr_bin(i,:) = mean(base_fr);
    base_mean(i,1) = mean(base_fr_bin(i,:));
    base_std(i,1) = std(mean(base_fr));
    thre(i,1) = base_mean(i,1) + 2*base_std(i,1);
    z_score_base(i,:) = (base_fr_bin(i,:) - base_mean(i,1))/base_std(i,1);
    
    for h = 1:size(range2,1)
        for k = 1:size(range2,2)-1
            beh_spike(h,k) = length(iTime(iTime >= range2(h,k) & iTime <= range2(h,k+1)));
            beh_fr(h,k) = beh_spike(h,k)/bins;
        end
    end
    mean_beh_fr(i,:) = mean(beh_fr);
    z_score_tmp = (mean_beh_fr(i,:) - base_mean(i,1))/base_std(i,1);
    z_score_judge = z_score_tmp(1:beh_judge/bins);
    z_score(i,:) = z_score_tmp;
    indx_up = find(z_score_judge > 2);
    diff_up = diff(indx_up);
    judge_up = isempty(find(diff_up == 1, 1));
    indx_down = find(z_score_judge < -2);
    diff_down = diff(indx_down);
    judge_down = isempty(find(diff_down == 1, 1));
   
    if judge_up == 0 && judge_down == 1
        z_score_pro(i,1) = {'excited'};
    else if judge_up == 1 && judge_down == 0
            z_score_pro(i,1) = {'inhibited'};
        else if judge_up == 0 && judge_down == 0
                error('max > 2 & min < -2')
            else z_score_pro(i,1) = {'unchanged'};
            end
        end
    end
end

    
cd(MatPath)
res.base_fr_bin = base_fr_bin;
res.z_score_base = z_score_base;
res.mean_beh_fr = mean_beh_fr;
res.z_score = z_score;
save([f1(1:end-4),'_fr.mat'],'res');

output = [base_mean,base_std,thre];
f2_new = reshape(f2,[],1);
chan_name = cellfun(@(x) x(1:end-4),f2_new,'UniformOutput',false);
title = {'filename','chanel', 'z_score_pro','base_mean','base_std','threshold','z_score'};
cd(ExcelPath)
xlswrite([f1(1:end-4),'.xlsx'],title,1,'A1');
filename = repmat({f1(4:end-4)},length(f2),1);
xlswrite([f1(1:end-4),'.xlsx'],filename,1,'A2'); 
xlswrite([f1(1:end-4),'.xlsx'],chan_name,1,'B2');
xlswrite([f1(1:end-4),'.xlsx'],z_score_pro,1,'C2');
xlswrite([f1(1:end-4),'.xlsx'],output,1,'D2');
xlswrite([f1(1:end-4),'.xlsx'],z_score,1,'G2');
close all
clear all
disp('finish')
%% part8 heatmap of ranking unit 
Matpath = ''; % fr_mat
bins = 0.25; 
bef_beh = 5;
beh_dur = 5;

filelist = sortdir(Matpath,'*home*.mat');
filelist2 = sortdir(Matpath,'*of*.mat');
cd(Matpath)
mycolor = load('ht_colormap.mat');
foxurine_zscore = [];
object_zscore = [];
for i = 1:length(filelist)
    filename = filelist{i};
    urine_tmp = load(filename);
    urine_matrix_tmp = [urine_tmp.res.z_score_base,urine_tmp.res.z_score];
    foxurine_zscore = [foxurine_zscore;urine_matrix_tmp];
end
foxurine_zscore_max = max(foxurine_zscore(:,21:28),[],2);

foxurine_new = [foxurine_zscore,foxurine_zscore_max];
foxurine_order_tmp = sortrows(foxurine_new,-size(foxurine_new,2));
foxurine_order_zscore = foxurine_order_tmp(:,1:end-1);
figure
imagesc(foxurine_order_zscore)
colormap(mycolor.CustomColormap);
caxis([-2 2]);
colorbar('location','East')
time_xaxis1 = 0.5:(1/bins):size(foxurine_order_zscore,2)+0.5;
time_xaxis1_label = -bef_beh:beh_dur;
set(gca,'XTick',time_xaxis1)
set(gca,'XTickLabel',time_xaxis1_label);
hold on
plot([20.5,20.5],[0.5 (size(foxurine_order_zscore,1)+0.5)],'k','linewidth',2);
set(gcf, 'position', get(0,'ScreenSize'));
title('urine','fontsize',15)
saveas(gcf,'urine_summary','bmp'); 
print(gcf,'-dpdf','-r600','urine_summary.pdf')



for i = 1:length(filelist2)
    filename2 = filelist2{i};
    object_tmp = load(filename2);
    object_matrix_tmp = [object_tmp.res.z_score_base,object_tmp.res.z_score];
    object_zscore = [object_zscore;object_matrix_tmp];
end
object_zscore_max = max(object_zscore(:,21:28),[],2);

object_new = [object_zscore,object_zscore_max];
object_order_tmp = sortrows(object_new,-size(object_new,2));
object_order_zscore = object_order_tmp(:,1:end-1);
figure
imagesc(object_order_zscore)
colormap(mycolor.CustomColormap);
caxis([-2 2]);
colorbar('location','East')
time_xaxis2 = 0.5:(1/bins):size(object_order_zscore,2)+0.5;
time_xaxis2_label = -bef_beh:beh_dur;
set(gca,'XTick',time_xaxis2)
set(gca,'XTickLabel',time_xaxis2_label);
hold on
plot([20.5,20.5],[0.5 (size(object_order_zscore,1)+0.5)],'k','linewidth',2);
title('object','fontsize',15)
saveas(gcf,'object_summary','bmp'); 
print(gcf,'-dpdf','-r600','object_summary.pdf')
clear all
close all
disp('finish')
    
    
    
    
    
    
    