clc
close all
clear

%% seed 持失 - MCS
mkdir([pwd,'\s_save_200228'])

nexp=50;
for trial=1:nexp
%     rng('default')
    file_name=[pwd,'\s_save_200228','\',sprintf('s_save_trial%d_200228.mat',trial)];
    s=rng('shuffle');%
    save(file_name,'s')
end

%% seed 持失 - Exp. Data
mkdir([pwd,'\s_save_exp_200228'])
nexp=50;
for trial=1:nexp
%     rng('default')
    file_name=[pwd,'\s_save_exp_200228','\',sprintf('s_save_exp_trial%d_200228.mat',trial)];
    s=rng('shuffle');%
    save(file_name,'s')
end

%% seed 持失 - Exp. Data for Sub Model
mkdir([pwd,'\s_save_sub_model_200228'])
nexp=50;
for trial=1:nexp
%     rng('default')
    file_name=[pwd,'\s_save_sub_model_200228','\',sprintf('s_save_exp_trial%d_200228.mat',trial)];
    s=rng('shuffle');%
    save(file_name,'s')
end