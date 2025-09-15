function [man_adv_x, man_adv_y] = CELLindivmotion1(statfile,thresh);




% derived from my CACTI-CI code (.../MATLAB/CELLTRACK_statdata.m)

% this code reads in cell stats data to track cell cores, then calculates
% an unique individual mean cell motion for each cell.

% cellfile is the stats file
% thresh = absolute val of the max u,v motion (m/s) you'll allow - designed to nan-out 'jumps'


% statfile =  cellstats
% thresh = 40;

%ncdisp(statfile)

trackbasetime = ncread(statfile,'base_time');
core_mean_y = single( ncread(statfile,'core_mean_y') );
core_mean_x = single( ncread(statfile,'core_mean_x') );

%  track_duration = single(ncread(statfile,'track_duration') );


[aa bb] = size(core_mean_x);

%%%% Manually calculate advection velocity from cell (x,y) as an alternative to teh mean adv_x,adv_y per radar scene used in FLEXTRKR to match cells between times:

man_adv_x = single(core_mean_y);   man_adv_x(:) = NaN; man_adv_y = man_adv_x;

for n = 1:bb
    for t = 2:aa
        if( isnan(core_mean_x(t,n)) == 0 & isnan(core_mean_x(t-1,n)) == 0)
            
            man_adv_x(t,n) = 1000*(  core_mean_x(t,n) - core_mean_x(t-1,n)  )/(trackbasetime(t,n)-trackbasetime(t-1,n));  %m/s
            man_adv_y(t,n) = 1000*(  core_mean_y(t,n) - core_mean_y(t-1,n)  )/(trackbasetime(t,n)-trackbasetime(t-1,n));  %m/s
    
        end

    end

end

%threshold out 'jumps'. I set this to 40 m/s for now, but 
%thresh = 40.0;  %m/s
man_adv_x(find(abs(man_adv_x) > thresh)) = NaN;
man_adv_y(find(abs(man_adv_y) > thresh)) = NaN;






% % smooth along time with window-size 2
% man_adv_x = movmean(man_adv_x,2,1,'omitnan');
% man_adv_y = movmean(man_adv_y,2,1,'omitnan');
% % kill the last non-nan time because it is artifical from the movmean:
% % note, this is hardcoded to a 2-sized window. You'll have to modify if you
% % change. 
% for n = 1:bb
%     nonnan = find(isnan(man_adv_x(:,n))==0);
%     if( isempty(nonnan)==0 )
%         man_adv_x(nonnan(end)) = NaN;
%         man_adv_y(nonnan(end)) = NaN;
%     end
% end
% 
% %set the frist time to the second time.
% man_adv_x(1,:) = man_adv_x(2,:);
% man_adv_y(1,:) = man_adv_y(2,:);



%return mean cell motion of first 45 min of cell
man_adv_x = mean(man_adv_x(1:3,:),1, 'omitnan');
man_adv_y = mean(man_adv_y(1:3,:),1, 'omitnan');









