function fscope_epm_zone_signal(MatPath,MatPath2,MatPath3,LogPath,ExcelPath,xlsname,framerate,filename)

cd(LogPath);
log = load(['log',filename,'.mat']);
signal_start_frame = log.markslog.framenum(2);
start_frame = log.markslog.framenum(3); 
finish_frame = log.markslog.framenum(4);

cd(MatPath)
info = load([filename,'.mat']);
signal = info.res.fscope_signal;
fscope_fps = info.res.fscope_fps;

signal_indx = 0:fscope_fps/framerate:length(signal);
for i = 1:(length(signal_indx)-1);
    new_signal(i,1) = mean(signal((signal_indx(i)+1):(signal_indx(i+1))));
end
beh_signal = new_signal((start_frame-signal_start_frame+1):(finish_frame-signal_start_frame+1));
baseline = mean(beh_signal); 
beh_signal_delF = (beh_signal - baseline)/baseline*100;

cd(MatPath2);
tmp1 = load([filename,'_mouse_point.mat']);
centr = tmp1.res.centr;
beh_centr = centr((start_frame:finish_frame),:);

tmp2 = load([filename,'_square.mat']);
ver_x = tmp2.ver_x;
ver_y = tmp2.ver_y;
hor_x = tmp2.hor_x;
hor_y = tmp2.hor_y;


%vertical
close1_ymin = 0;
close1_ymax = ver_y(1);
close2_ymin = ver_y(2);
close2_ymax = 600;
v_xmin = hor_x(1);
v_xmax = hor_x(2);

% horizontal
open1_xmin = 0;
open1_xmax = hor_x(1);
open2_xmin = hor_x(2);
open2_xmax = 600;
h_ymin = ver_y(1);
h_ymax = ver_y(2);

posi_indx = cell(1,5);

close1_up_indx = find(v_xmin <= beh_centr(:,1) & beh_centr(:,1) <= v_xmax & close1_ymin <= beh_centr(:,2) & beh_centr(:,2) < close1_ymax);
close2_down_indx = find(v_xmin <= beh_centr(:,1) & beh_centr(:,1) <= v_xmax & close2_ymin < beh_centr(:,2) & beh_centr(:,2) <= close2_ymax);
open1_left_indx = find(open1_xmin <= beh_centr(:,1) & beh_centr(:,1) < open1_xmax & h_ymin <= beh_centr(:,2) & beh_centr(:,2) <= h_ymax);
open2_right_indx = find(open2_xmin < beh_centr(:,1) & beh_centr(:,1) <= open2_xmax & h_ymin <= beh_centr(:,2) & beh_centr(:,2) <= h_ymax);
center_indx = (1:size(beh_centr,1))'; 
center_indx([close1_up_indx;close2_down_indx;open1_left_indx;open2_right_indx]) = [];

posi_indx{1,1} = close1_up_indx;
posi_indx{1,2} = close2_down_indx;
posi_indx{1,3} = open1_left_indx;
posi_indx{1,4} = open2_right_indx;
posi_indx{1,5} = center_indx;

for i = 1:size(posi_indx,2)
    if isempty(posi_indx{1,i}) == 1
        zone_count{1,i} = 0;
        zone_time{1,i} = 0;
        zone_meandur{1,i} = 0;
    else data_tmp = posi_indx{1,i};
        if length(data_tmp) == 1;
            zone_count{1,i} = 1;
            zone_time{1,i} = 1/framerate;
            zone_meandur{1,i} = 1/framerate;
            transfer_point{1,i} = data_tmp;
        else data_dif = diff(data_tmp);
             indx_tmp = find(data_dif > 25); 
             if isempty(indx_tmp) == 1;
                 zone_count{1,i} = 1;
                 zone_time{1,i} = length(data_tmp)/framerate;
                 zone_meandur{1,i} = length(data_tmp)/framerate;
                 transfer_point{1,i} = sort([data_tmp(1);data_tmp(end)]);
             else if isempty(indx_tmp) ~= 1 && indx_tmp(end)+1 == length(data_tmp);
                     indx_tmp(end) = [];
                     zone_count{1,i} = length(indx_tmp)+1;
                     zone_time{1,i} = (length(data_tmp)-1)/framerate;
                     zone_meandur{1,i} = zone_time{1,i}/zone_count{1,i};
                     indx_tmp2 = indx_tmp + 1;
                     transfer_point{1,i} = sort([data_tmp(1);data_tmp(indx_tmp);data_tmp(indx_tmp2);data_tmp(end-1)]);
                 else
                  zone_count{1,i} = length(indx_tmp)+1;
                  zone_time{1,i} = length(data_tmp)/framerate;
                  zone_meandur{1,i} = zone_time{1,i}/zone_count{1,i};
                  indx_tmp2 = indx_tmp + 1;
                  transfer_point{1,i} = sort([data_tmp(1);data_tmp(indx_tmp);data_tmp(indx_tmp2);data_tmp(end)]);
                 end
             end
        end
    end
    indx_tmp = [];
    indx_tmp2 = [];
    data_tmp = [];
end

for i = 1:size(posi_indx,2)
    beh_zone_signal{1,i} = beh_signal_delF(posi_indx{1,i});
    zone_signal{1,i} = mean(beh_zone_signal{1,i});
end

close_time = zone_time{1,1} + zone_time{1,2};
close_entry = zone_count{1,1} + zone_count{1,2};
close_meandur = close_time/close_entry;
open_time = zone_time{1,3} + zone_time{1,4};
open_entry = zone_count{1,3} + zone_count{1,4};
open_meandur = open_time/open_entry;

close_signal_point = [beh_zone_signal{1,1};beh_zone_signal{1,2}];
open_signal_point = [beh_zone_signal{1,3};beh_zone_signal{1,4}];
center_signal_point = beh_zone_signal{1,5};
close_signal = mean(close_signal_point);
open_signal = mean(open_signal_point);
center_signal = mean(center_signal_point);

cd(ExcelPath);
title = {'filename','mouse','virus','behavior','trial','close1 time','close2 time','open1 time','open2 time','center time','close1 count',...
    'close2 count','open1 count','open2 count','center count','close1 meandur','close2 meandur','open1 meandur','open2 meandur','center meandur',...
    'close1 signal','close2 signal','open1 signal','open2 signal','center signal','close time','close count','close menadur','open time','open entry',...
    'open menadur','close signal','open signal'};
xlswrite(xlsname,title,1,'A1');

indx = strfind(filename,'_');
excel{1,1} = filename;
excel{1,2} = filename(2:indx(1)-1);
excel{1,3} = filename(indx(2)+1:indx(3)-1);
excel{1,4} = filename(indx(4)+1:indx(5)-1);
excel{1,5} = filename(indx(5)+1:indx(6)-1);
excel{1,6} = zone_time{1};
excel{1,7} = zone_time{2};
excel{1,8} = zone_time{3};
excel{1,9} = zone_time{4};
excel{1,10} = zone_time{5};
excel{1,11} = zone_count{1};
excel{1,12} = zone_count{2};
excel{1,13} = zone_count{3};
excel{1,14} = zone_count{4};
excel{1,15} = zone_count{5};
excel{1,16} = zone_meandur{1};
excel{1,17} = zone_meandur{2};
excel{1,18} = zone_meandur{3};
excel{1,19} = zone_meandur{4};
excel{1,20} = zone_meandur{5};
excel{1,21} = zone_signal{1};
excel{1,22} = zone_signal{2};
excel{1,23} = zone_signal{3};
excel{1,24} = zone_signal{4};
excel{1,25} = zone_signal{5};
excel{1,26} = close_time;
excel{1,27} = close_entry;
excel{1,28} = close_meandur;
excel{1,29} = open_time;
excel{1,30} = open_entry;
excel{1,31} = open_meandur;
excel{1,32} = close_signal;
excel{1,33} = open_signal;

[num tmp] = xlsread(xlsname,1, 'A1:A10000');
xlswrite(xlsname,excel,1,['A',num2str(size(tmp,1)+1)]);

for i = 1:size(transfer_point,2)
    if isempty(transfer_point{1,i}) == 1
        i = i+1;
    else transfer_point_tmp = transfer_point{1,i};
        transfer_in{1,i} = transfer_point_tmp(1:2:length(transfer_point_tmp));
        transfer_out{1,i} = transfer_point_tmp(2:2:length(transfer_point_tmp));
    end
end

close_transin_frame = sort([transfer_in{1};transfer_in{2}]);
close_transout_frame = sort([transfer_out{1};transfer_out{2}]);
open_transin_frame = sort([transfer_in{3};transfer_in{4}]);
open_transout_frame = sort([transfer_out{3};transfer_out{4}]);
center_transin_frame = transfer_in{5};
center_transout_frame = transfer_out{5};

cd(MatPath3)
res.beh_signal = beh_signal;
res.beh_signal_delF = beh_signal_delF;
res.close_signal_point = close_signal_point;
res.open_signal_point = open_signal_point;
res.center_signal_point = center_signal_point;
res.close_transin_frame = close_transin_frame;
res.close_transout_frame = close_transout_frame;
res.open_transin_frame = open_transin_frame;
res.open_transout_frame = open_transout_frame;
res.center_transin_frame = center_transin_frame;
res.center_transout_frame = center_transout_frame;    
save([filename,'_epm_beh_info'],'res');
clear all
close all