function fscope_align_sigal_EPM_heatmap(FscopeMatPath,MatPath,LogPath,FigPath,FileName,nbins,beh_time,color_caxis)


%% load the start point of the signal recording and the behavior start
cd(LogPath)
scorelog = load(['log',FileName,'.mat']);
marksName = scorelog.markslog.framename;
marksFrame = scorelog.markslog.framenum;
if size(marksName,1)>3
    signal_start_frame = marksFrame(2);
    beh_start_frame = marksFrame(3);
end
%% load the track of the mouse 
cd(MatPath)
load([FileName,'_mouse_point.mat']);
centr = res.centr;
width = res.width;
height = res.height;
FrameRate = res.FrameRate;
Img11 = [];
clear track;
centr = round(centr((beh_start_frame:(beh_start_frame + beh_time*60*FrameRate)),:));
for i = 1:width/nbins
    for j = 1:height/nbins
        idx = find(1+(i-1)*nbins <= centr(:,1) & centr(:,1) <= nbins*i & 1+(j-1)*nbins <= centr(:,2) & centr(:,2) <= nbins*j);%找到在同一个格子里的帧数
        Img11(j,i) = length(idx)/FrameRate;
        data_indx = find_continuous(idx); 
        track{j,i} = data_indx;
    end
    
end
% Img_log = log(Img11);
% cd(FigPath)
% figure;imagesc(Img_log)
% colorbar
% % caxis(color_caxis);
% caxis(color_caxis);
% title([FileName,'_track_map_bins = ',num2str(nbins)],'FontSize',10,'Interpreter','none');
% set(gcf,'position', get(0,'ScreenSize'));
% saveas(gcf, [FileName,'_track_map_bins =',num2str(nbins),'.tif'],'tif');

%% load the fscope data % normalized by z-score deltaF/F
cd(FscopeMatPath)
load([FileName,'.mat']);
fscope_fps  = res.fscope_fps;
fscope_signal_all = res.fscope_signal;
behavior_start_sample = (beh_start_frame-signal_start_frame)/FrameRate*fscope_fps+1;
behavior_end_sample = (beh_start_frame+beh_time*60*FrameRate-signal_start_frame)/FrameRate*fscope_fps+1; 
fscope_signal = fscope_signal_all(behavior_start_sample:behavior_end_sample);
mean_signal = mean(fscope_signal);
fscope_signal = (fscope_signal - mean_signal)/mean_signal*100; 
clear ave_tot;
clear max_tot;
clear median_tot;
for i = 1:size(track,1)
    for j = 1:size(track,2)
        tmp = track{i,j};
        tmp_point = round((tmp-1)/FrameRate*fscope_fps+1);
        ave_signal = zeros(size(tmp_point,1),1);
        for hh = 1:size(tmp_point,1)
        signal_tmp = fscope_signal(tmp_point(hh,1):tmp_point(hh,2));
        ave_signal(hh) = mean(signal_tmp);
        end
        ave_tot(i,j) = mean(ave_signal);
    end
end
cd(FigPath)
figure;
imagesc(ave_tot);colorbar;
caxis(color_caxis);
title([FileName,'_deltaF_F_ave_bins = ',num2str(nbins)],'FontSize',10,'Interpreter','none');
set(gcf,'position', get(0,'ScreenSize'));
saveas(gcf, [FileName,'_deltaF_F_ave_bins = ',num2str(nbins)],'tif');
print(gcf,'-dpdf','-r600',[FileName,'_deltaF_F_ave_bins = ',num2str(nbins),'.pdf']);

close all;

%%
res.signal_start_frame = signal_start_frame;
res.beh_start_frame = beh_start_frame;
res.nbins = nbins;
res.Img_log = Img_log;
res.beh_time = beh_time;
res.ave_zscore_tot = ave_tot;
res.track = track; %
res.fscope_signal = fscope_signal;%
cd(MatPath)
save([FileName,'_deltaF_F_signal_nbins = ',num2str(nbins),'.mat'],'res');
