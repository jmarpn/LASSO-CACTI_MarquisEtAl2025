function [layermean]  = LAYERMEAN(var,z,kbot,ktop);

% %diagnostics:
% kbot = 189;
% ktop = 210;
% z = ht;
%var = rh(:,1000);

%var(kbot) = 100
%  var(kbot:ktop)
%  z(kbot:ktop)

%var = u_crel;

%chop out NaNs:
keepervar = var(kbot:ktop);
keeperz = z(kbot:ktop);
kill = find(isnan(keepervar));
keepervar(kill) = [];
keeperz(kill) = [];

layermean = 0;
for k = 1:length(keepervar)-1
    clear chunk avg dz
    avg = mean(keepervar(k:k+1),'omitnan')          ;
    dz = (keeperz(k+1) - keeperz(k))                      ; 
    chunk = avg * dz./( keeperz(end) - keeperz(1) )        ;
    
    layermean = layermean + chunk             ; 
end

if( length(find(isnan(keepervar))) == length(keepervar) | isempty(keepervar) == 1 )
    layermean(:) = NaN;
end




% layermean = 0;
% for k = kbot:ktop-1
%     clear chunk avg dz
%     avg = mean(var(k:k+1),'omitnan')          ;
%     dz = (z(k+1) - z(k))                      ; 
%     chunk = avg * dz./( z(ktop) - z(kbot) )        ;
%     
%     layermean = layermean + chunk             ; 
% end
%layermean

%mean(var(kbot:ktop))