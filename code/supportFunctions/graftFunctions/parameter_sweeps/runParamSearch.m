poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end
if ~ispc
    parpool(12)
else
    parpool(2)
end

%%
jobs = importParams('jobs.txt');
jobs = [jobs, zeros(size(jobs,1), 3)];
parfor i = 1:size(jobs,1)
    fprintf('%d \n',i);
    tStart = tic;
    job = jobs(i,:);
    lambda1    = job(3);
    lambda2    = job(4);
    lambda3    = job(5);
    hiddenSize = job(1);
    nEigs      = job(2);
    %%
    
    [acc_train, acc_test] = runScript(lambda1,lambda2,lambda3,hiddenSize,nEigs);
    
    job(6)  = acc_train;
    job(7)  = acc_test;
    job(8) = toc(tStart);
    
    fid=fopen('out.txt','a');
    fprintf(fid,'param %d \t %d \t %d \t %6.2e \t %6.2e \t %6.2e \t %6.4f \t %6.4f \t %3.2f \r\n', ...
        i,job(1),job(2),job(3),job(4),job(5),job(6),job(7),job(8));
    fclose(fid);
    jobs(i,:) = job;
    
end

save('output.mat','jobs');
