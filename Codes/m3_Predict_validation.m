%% cross dimension validation
nsub = 34;
nlevel = 4;
G1_data = load('ART_G1_4lvl_predictiondata.mat');
G2_data = load('ART_G2_4lvl_predictiondata.mat');

%% cross prediction based on 10-fold cross-validation
acc_inter=zeros(2,10,10);
acc_within=zeros(2,10);
for repeat = 1:10
    CVindex = GenerateCV(nsub, nlevel, repeat,10);
    [~, stats_G1] = predict(G1_data.AATpredition,  'algorithm_name', 'cv_svr', 'nfolds', CVindex, 'error_type', 'mse');
    [~, stats_G2] = predict(G2_data.AATpredition,  'algorithm_name', 'cv_svr', 'nfolds', CVindex, 'error_type', 'mse');
    acc_within(1,repeat) = stats_G1.pred_outcome_r;
    acc_within(2,repeat) = stats_G2.pred_outcome_r;
    for cv_i = 1:10
        cv_Y_a = G1_data.AATpredition.Y(CVindex==cv_i);
        dat_subj_a = G1_data.AATpredition.dat(:, CVindex==cv_i);
        cv_Y_c = G2_data.AATpredition.Y(CVindex==cv_i);
        dat_subj_c = G2_data.AATpredition.dat(:, CVindex==cv_i);
        
        a_weights =stats_G1.other_output_cv{cv_i,1}(:,1);
        cv_a_c = a_weights' * dat_subj_c;
        acc_inter(1,cv_i,repeat) = corr(cv_Y_c,cv_a_c');

        c_weights =stats_G2.other_output_cv{cv_i,1}(:,1);
        cv_c_a = c_weights' * dat_subj_a;
        acc_inter(2,cv_i,repeat) = corr(cv_Y_a,cv_c_a');
    end
end

cross_repeats =squeeze(mean(acc_inter,2))';
in_repeats = acc_within';

crossprediction_data_plot = [cross_repeats(:,1),in_repeats,cross_repeats(:,2)];
boxplot([cross_repeats(:,1),in_repeats,cross_repeats(:,2)])

% test the indivdiual predictive performance between models 
[h,p,ci,stats] =ttest(crossprediction_data_plot(:,1),crossprediction_data_plot(:,2))
[h,p,ci,stats] =ttest(crossprediction_data_plot(:,3),crossprediction_data_plot(:,4))