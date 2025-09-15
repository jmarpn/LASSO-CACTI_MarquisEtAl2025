

cluster = parcluster('local');
cluster.NumWorkers(128)
saveAsProfile(cluster,'perl_cluster')
pp = parpool('perl_cluster', 128)
  
delete(gcp('nocreate'))
  
pc = parcluster('local')
parpool(pc, 128);
  
spmd rank = labindex;
    fprintf(1,'Hello from %d\n',rank);
end

parfor i = 1:128 
    rank=i; 
    fprintf(1,'Hello from core %d\n',rank); 
end 

delete(gcp('nocreate'));





