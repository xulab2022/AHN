function fscope_of_no_combine_during_record_zone_signal(FscopeMatpath,MousepointMatpath,LogPath,ExcelPath,xlsname,filename,video_fps,zone,beh_time)
% FscopeMatpath = '';
% MousepointMatpath = '';
% LogPath = ''; 
% ExcelPath = '';
% xlsname = '';
% filename = '';
% video_fps = 25; %25frames/s
% zone = [20,30,40]; %cm
% beh_time = 600;%%% s

cd(LogPath);
info1 = load(['log',filename,'.mat']);
signal_start_frame = info1.markslog.framenum(2);
beh_start_frame = info1.markslog.framenum(3);
beh_end_frame = beh_start_frame+beh_time*video_fps;


cd(FscopeMatpath);
info2 = load([filename,'.mat']);
fscope_fps = info2.res.fscope_fps;
signal = info2.res.fscope_signal;

signal_indx = 0:fscope_fps/video_fps:length(signal);
for i = 1:(length(signal_indx)-1);
    new_signal(i,1) = mean(signal((signal_indx(i)+1):(signal_indx(i+1)))); %signal transfer to value of per frame
end


beh_sig{1,1} = new_signal(1:1+beh_time*video_fps);
beh_sig{1,2} = new_signal(beh_start_frame-signal_start_frame+1:beh_end_frame-signal_start_frame+1);
baseline = mean(beh_sig{1,1});
baseline2 = mean(beh_sig{1,2});
for i = 1:size(beh_sig,2)
    delta_signal_beh{1,i} = (beh_sig{1,i} - baseline)/baseline*100; 
end


cd(MousepointMatpath);
info3 = load([filename,'_mouse_point.mat']);
raw_centr = info3.res.centr;
beh_centr{1,1} = raw_centr((signal_start_frame : signal_start_frame+beh_time*video_fps),:);
beh_centr{1,2} = raw_centr((beh_start_frame : beh_end_frame),:);


zone_min = round(300*(1-zone/40));
zone_max = round(300*(1+zone/40));

for i = 1:size(beh_centr,2)
    beh_centr_tmp = beh_centr{1,i};
    for j = 1:length(zone)
        position_indx_tmp{1,j} = find(zone_min(j)<=beh_centr_tmp(:,1) & beh_centr_tmp(:,1)<= zone_max(j) & zone_min(j)<=beh_centr_tmp(:,2) & beh_centr_tmp(:,2)<=zone_max(j));
    end
    position_indx = position_indx_tmp;
    for z = 1:size(position_indx_tmp,2) - 1
        same_indx_tmp = ismember(cell2mat(position_indx_tmp(z+1)),cell2mat(position_indx_tmp(z)));
        same_indx = same_indx_tmp == 1;
        position_indx_tmp_tmp = cell2mat(position_indx_tmp(z+1));
        position_indx_tmp_tmp(same_indx) = [];
        position_indx{z+1} = position_indx_tmp_tmp;
        same_indx_tmp = [];
        same_indx = [];
        position_indx_tmp_tmp = [];
    end
    beh_position_indx{1,i} =  position_indx;
    position_indx_tmp = [];
    position_indx = [];
    beh_centr_tmp = [];
end
    
for i = 1:size(beh_position_indx,2);
    beh_position_indx_tmp = beh_position_indx{1,i};
    for j = 1:size(beh_position_indx_tmp,2)
        if isempty(beh_position_indx_tmp{1,j}) == 1;
            mean_signal(1,j) = NaN;
            max_signal(1,j) = NaN;
        else
            position_signal{1,j} = delta_signal_beh{1,i}(beh_position_indx_tmp{1,j});
            mean_signal(1,j) = mean(position_signal{1,j});
            max_signal(1,j) = max(position_signal{1,j});
            time(1,j) = length(position_signal{1,j})/video_fps;
        end
    end
    zone_mean_sig{1,i} = mean_signal;
    zone_max_sig{1,i} = max_signal;
    zone_time{1,i} = time;
    mean_signal = [];
    max_signal = [];
    time = [];
end

for i = 1:size(beh_position_indx,2)
    beh_position_indx_tmp2 = beh_position_indx{1,i};
    for j = 1:size(beh_position_indx_tmp2,2)
        position_frame = beh_position_indx_tmp2{1,j};
        position_intval = diff(position_frame);
        position_number = find(position_intval ~= 1);
        position_entry(1,j) = length(position_number) + 1;
        position_frame = [];
        position_intval = [];
        position_number = [];
    end
    zone_entry{1,i} = position_entry;
end
 

    
cd(ExcelPath)
title1 = {'filenmae','mouse','virus','behavior','trial','of baseline','no baseline','of_center_signal','of_middle_signal','of_periphery_signal','no_center_signal','no_middle_signal','no_periphery_signal','of_center_time','of_middle_time','of_periphery_time',...
    'no_center_time','no_middle_time','no_periphery_time','of_max_center_signal','of_max_middle_signal','of_max_periphery_signal','no_max_center_signal','no_max_middle_signal','no_max_perophery_signal','of_center_entry','of_middle_entry','of_periphery_entry',...
    'no_center_entry','no_middle_entry','no_periphery_entry','of_center_mean_dur','of_middle_mean_dur','of_periphery_mean_dur','no_center_mean_dur','no_middle_mean_dur','no_periphery_mean_dur'}; %%%%% for of use 3 part
xlswrite(xlsname,title1,1,'A1');
indx = strfind(filename,'_');
excel{1,1} = filename;
excel{1,2} = filename(2:indx(1)-1);
excel{1,3} = filename(indx(2)+1:indx(3)-1);
excel{1,4} = filename(indx(4)+1:indx(5)-1);
excel{1,5} = filename(indx(5)+1:indx(6)-1);
excel{1,6} = baseline;
excel{1,7} = baseline2;

tmp1 = cell2mat(zone_mean_sig);
tmp2 = cell2mat(zone_time);
tmp3 = cell2mat(zone_max_sig);
tmp4 = cell2mat(zone_entry);
tmp5 = tmp2./tmp4;

for i = 1:length(tmp1)
    excel{1,7+i} = tmp1(i);
    excel{1,7+length(tmp1)+i} = tmp2(i);
    excel{1,7+2*length(tmp1)+i} = tmp3(i);
    excel{1,7+3*length(tmp1)+i} = tmp4(i);
    excel{1,7+4*length(tmp1)+i} = tmp5(i);
end

[num,tmp] = xlsread(xlsname,1);
xlswrite(xlsname,excel,1,['A',num2str(size(tmp,1)+1)]);
clear all
close all