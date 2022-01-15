function makeParamSweepData(param_set_dir,param_set_file)

% makeParamSweepData(param_set_dir,param_set_file)
% 
% 
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create data

N_reps                      = 10;
sim_sweep_opts.noise_var    = linspace(0,0.1,5);
sim_sweep_opts.nBG          = [0,1];

dat_fields = fieldnames(sim_sweep_opts);
data_szs   = zeros(numel(dat_fields),1);
for kk = 1:numel(dat_fields)
    data_szs(kk) = numel(sim_sweep_opts.(dat_fields{kk}));
end
N_datasets = prod(data_szs)*N_reps;                                        % Record the number of datasets that will be generated

sim_opts.dims      = [50,50,400];
sim_opts.n_dict    = 13;
sim_opts.poles     = -0.7;
sim_opts.p_evt     = 0.95;
sim_opts.noise_var = 0.001;
sim_opts.AR1       = false;
sim_opts.prof_type = 'blob';
sim_opts.nBG       = 1;

Ix = cell(size(data_szs));                                                 % Initialize an array to store the data parameter set to use here
for kk = 1:prod(data_szs)
    [Ix{:}] = ind2sub(data_szs, kk);
    for ll = 1:numel(dat_fields)
        sim_opts.(dat_fields{ll}) = sim_sweep_opts.(dat_fields{ll})(Ix{ll});
    end
    for ll = 1:N_reps
        data_num       = kk*(N_reps-1)+ll;
        [sim_cube,P,D] = simSpatialData(sim_opts);
        save([param_set_dir,param_set_file,'_sim',sprintf('%05.0f',data_num),'.mat'],'sim_cube','P','D')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up param struct

param_struct.lambda    = linspace(0.1,10,3);
param_struct.lambda2   = linspace(0.1,0.5,1);
param_struct.lambda4   = linspace(0.1,0.5,1);
param_struct.beta      = [0.001,0.01,0.1];
param_struct.n_dict    = [10,13,20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up defaults for the other variables


params.lambda3   = 0.5;
%params.n_dict    = 20;
params.tau       = 1;
params.learn_eps = 1e-4;
params.max_learn = 50;
params.grad_type = 'full_ls_cor';
params.nonneg    = true; 
params.nneg_dict = true; 

corr_kern.w_time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the parameters

par_fields = fieldnames(param_struct);
par_szs    = zeros(numel(par_fields),1);
for kk = 1:numel(par_fields)
    par_szs(kk) = numel(param_struct.(par_fields{kk}));
end
N_params = prod(par_szs);
N_jobs   = N_params*N_datasets;
save([param_set_dir,param_set_file,'.mat'],'param_struct','params','corr_kern','N_datasets','N_params','N_jobs')

fprintf('Created dataset with %d datasets, %d parameter sets for a total of %d jobs.',N_datasets,N_params,N_jobs)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
