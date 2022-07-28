function [fname ftime] = sortdir(path2sort, nameclue, bytime)
%SORTDIR Sorted file list of a dir, can also sort by
%    time. 
%   [fname ftime] = SORTDIR(path2sort, nameclue, [bytime]) 
%   
%   Input:
%       PATH2SORT(str) - 
%       NAMECLUE(str) - e.x. '*.mat' 
%       BYTIME(bool) - true is sort by time
%   Output:
%       FNAME(cell) - 
%       FTIME(cell) - 
% 
%   Notes: think about whether to use lower before sort
% 
%   See also 
% 
%   by Cheng Wang (chengwang@ion.ac.cn), 2010-06-01.

if ~exist('bytime', 'var')
    bytime = false;
end

tfname = dir(fullfile(path2sort, nameclue));

if ~bytime
    fname = sortrows({tfname.name}');
else
    ftime = cell2mat(cellfun(@datevec, {tfname.date}', 'UniformOutput', ...
                    false));
    fname = {tfname.name}';
    [ftime,idx] = sortrows(ftime);
    fname = fname(idx);
end

return;
