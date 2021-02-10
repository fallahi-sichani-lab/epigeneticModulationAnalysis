
%% epigeneticModulationAnalysis
% Epigenetic modulation reveals differentiation state specificity of oncogene addiction

% Custom MATLAB scripts for the analysis of drug response data, as presented in the manuscript:
% Khaliq M, Manikkam M, Martinez ED, and Fallahi-Sichani M. Epigenetic modulation reveals differentiation state specificity of oncogene addiction.

%% Experimental conditions

cell_lines = [
    "Primary melanocytes"
    "MMACSF"
    "C32"
    "WM902B"
    "WM115"
    "WM2664"
    "UACC62"
    "SKMEL19"
    "SKMEL28"
    "WM1552C"
    "LOXIMVI"
    "HS294T"
    "RVH421"
    "COLO858"
    "A2058"
    "A375"
    "A375-NRAS(Q61K)"
    ];

epigenetic_compounds = [
    "Control"
    "Givinostat"
    "CUDC-907"
    "JIB-04"
    "AZ6102"
    "I-BET-762"
    "OTX015"
    "SP2509"
    ];

timepoints = [
    0
    72
    120
    ];

MAPK_drugs = [
    "DMSO"
    "Vemurafenib (100 nM)"
    "Vemurafenib (100 nM) + Trametinib (10 nM)"
    ];
    

%% Compute "net growth rate" from log_normalized_cell_number, i.e. cell count data normalized to the first timepoint

clear net_growth_rate time_avg_net_growth_rate;

for cell_line_id = 1:length(cell_lines)
    for epigenetic_compound_id = 1:length(epigenetic_compounds)
        for epigenetic_dose_id = 1:2
            for MAPK_drug_id = 1:length(MAPK_drugs)
                for timepoint_id = 1:length(timepoints)
                    if (timepoint_id == 1)
                        net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id,timepoint_id) = nan;
                    else
                        
                        net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id,timepoint_id) = ...
                            (mean_log_normalized_cell_number(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id,timepoint_id)-...
                            mean_log_normalized_cell_number(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,1,1)) / ...
                            ((timepoints(timepoint_id)-timepoints(1))/24);
                    end
                end
                % time-averaged growth rate
                time_avg_net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                    nanmean(net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id,:));
            end
        end
    end
end

%% Calculate combination effectiveness based on Bliss Independence

clear epi_drug_combined_effect; % This is to define the comined effect of drugs individually or in combination (based on DIP)
clear epi_drug_combined_effect_fa; % Fraction of cells affected based lon DIP
clear epi_drug_combined_effect_Bliss_metric; % This is defined to evaluate combination effectiveness using Bliss Independence (based on DIP metric)

% min_DIP represents the minimum of DIP values reported across all cell lines and drug treatment conditions in this study.
min_DIP = -4.5528;

for MAPK_drug_id = 1:length(MAPK_drugs)
    for cell_line_id = 1:length(cell_lines)
        for epigenetic_compound_id = 1:length(epigenetic_compounds)
            for epigenetic_dose_id = 1:2
                if MAPK_drug_id == 1
                    % definition of DIP (based on Harris at al.)
                    epi_drug_combined_effect(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                        (time_avg_net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id))...
                        / time_avg_net_growth_rate(cell_line_id,1,1,1);
                    epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                        (1 - epi_drug_combined_effect(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id)) / ...
                        (1-(min_DIP));
                    
                    % definition of Bliss Independence model based on DIP
                    epi_drug_combined_effect_Bliss_metric(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = nan;
                    
                else
                    % definition of DIP (based on Harris at al.)
                    epi_drug_combined_effect(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                        (time_avg_net_growth_rate(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id))...
                        / time_avg_net_growth_rate(cell_line_id,1,1,1);
                    epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                        (1 - epi_drug_combined_effect(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id)) / ...
                        (1-(min_DIP));
                    if epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) < 0
                        epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = 0;
                    elseif epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) > 1
                        epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = 1;
                    end
                    
                    
                    %definition of Bliss Independence model based on DIP
                    epi_drug_combined_effect_Bliss_metric(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id) = ...
                        (epi_drug_combined_effect_fa(cell_line_id,1,epigenetic_dose_id,MAPK_drug_id) + epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,epigenetic_dose_id,1) -...
                        epi_drug_combined_effect_fa(cell_line_id,1,epigenetic_dose_id,MAPK_drug_id) * epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,1)) /...
                        epi_drug_combined_effect_fa(cell_line_id,epigenetic_compound_id,epigenetic_dose_id,MAPK_drug_id);
                end
            end
        end
    end
end

%% Perform t-SNE analysus on log(Integrated_single_cell_data) (log-transformed multiplexed single-cell data, e.g., Mitf, Ngfr, Axl)

iteration_factor = 2000;
Perplexity_factor = 480;
options = statset('MaxIter',iteration_factor);
rng default;

[pca_coeff,pca_score,pca_latent,pca_tsquared,pca_explained,pca_mu] = pca(zscore(log(Integrated_single_cell_data)));
clear tsne_Integrated_single_cell_data;
tsne_Integrated_single_cell_data = tsne(pca_score(:,1:2),'Algorithm','barneshut','Perplexity',Perplexity_factor,'Options',options,'Exaggeration',4,'LearnRate',1000);


%% PLSR analysis using input matrix (X; z-scored matrix of signaling data) and ouput vector (Y; z-scored vector of drug response data)

clear TSS;
clear PLSR_XLoading PLSR_YLoading PLSR_XScore PLSR_YScore PLSR_yfit;
clear Rsquare Qsquare;

[n_observations n_variables] = size(X);
TSS = sum((Y-mean(Y)).^2);
ncomp = min(8,length(Y)-2);
clear XLoading YLoading XScore YScore BETA PCTVAR MSE stats;
[XLoading,YLoading,XScore,YScore,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp,'cv',length(Y),'mcreps',100);

% Prediction accuracy (leave-one-out cross validation)
Qsquare = [0 1-length(Y)*MSE(2,2:end)/TSS];
% Performance
Rsquare = [0 cumsum(PCTVAR(2,:))];

% Calculate VIP (Variable Importance in Projection) Scores
[MMM III] = max(squeeze(Qsquare(2:end)));
optimized_ncomp_for_VIP = min(III,3);

if Qsquare(optimized_ncomp_for_VIP+1) < 0
    optimized_ncomp_for_VIP = 0;
end


sum1 = zeros(1,n_variables);
sum2 = 0;
clear SS Wnorm2;
for i = 1:optimized_ncomp_for_VIP % number of principle components
    SS(i) = (YLoading(i)^2)*(XScore(:,i)'*XScore(:,i));
end
for i = 1:optimized_ncomp_for_VIP
    sum2 = sum2 + SS(i);
    Wnorm2(i) = stats.W(:,i)'*stats.W(:,i);
end

clear VIP;
for counter = 1:n_variables
    for k = 1:optimized_ncomp_for_VIP
        sum1(counter) = sum1(counter) + SS(k)*stats.W(counter,k)^2/Wnorm2(k);
    end
    VIP(counter) = (n_variables*sum1(counter)/sum2)^0.5;
end


%% MLR analysis using input matrix (X; z-scored matrix of signaling data) and ouput vector (Y; z-scored vector of drug response data)

[b,bint,r,rint,stats] = regress(Y,X);
Yfit = X*b;
for i = 1:length(Y)
    plot(Y,Yfit,'marker','o','markersize',10,'MarkerEdgeColor','k','linestyle','none');
end
xlabel('Measured norm growth rate');
ylabel('Predicted norm growth rate');
[pearson_corr pearson_pvalue] = corr(Y,Yfit,'rows','complete')
