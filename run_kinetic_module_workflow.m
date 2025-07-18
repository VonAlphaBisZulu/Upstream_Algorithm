function run_all_upstream_MM(f)
disp(f)
% add location of required folders and tools
addpath(genpath('Code/'))

%% ! this part has to be modified by user 
addpath(genpath('../Git/cobratoolbox')) % path to cobra toolbox
addpath(genpath('../F2C2')) % path to F2C2 
%pathToR = 'C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe'; % path to R, required packages igraph, R.matlab

% run code for ArabidopsisCoreModel as an example
files = dir('Models_original/*.xml'); % location of model(s), can be .sbml format as well
% index of the model file in files for which the code should be run 

%% start with preprocessing and spliting reactions to elementary reaction format
[model,name] = run_preprocessing(files,f);

%% find balanced complexes
disp('balanced')
[B,group] = find_balanced_complexes(model);

%% sampling to reduce the run time of finding concordant complexes

% we use complex x complex matrix CC to store the found relationships
% 1 balanced complex (diagonal), -1 trivially concordant, 2 - non-trivial concordant

% At - list of complex pairs to check for concordance
disp('sampling')
num_samples = 10; % number of samples per round 
[At,CC] = pre_concordant_sampling(model,name,B,num_samples);

%% find concordant complexes
% to allow that this part is run seperately for subsets of pairs to check
% we used start and stop valiables to define the range to check in one job
% later the CC matrices can then be combined
start = 1; stop = size(At,1);
disp('concordant')
CC = find_concordant_complexes(model,group,CC,At,start,stop);

Results_balanced.MODEL_r{1} = model; % renaming to run code for kinetic modules 
    
% Group mutually concordant complexes with balanced complexes

[CP(:,1),CP(:,2)] = find(CC~=0);
unclassified = unique(CP);

class_with_balanced=[];
while ~isempty(unclassified)
    i = unclassified(1);
    class_with_balanced{end+1} = i;
    [r,~]=find(ismember(CP,i));
    while length(unique(reshape(CP(r,:),[],1))) > length(class_with_balanced{end})
        class_with_balanced{end} = unique(reshape(CP(r,:),[],1));
        unclassified = setdiff(unclassified,class_with_balanced{end});
        i = class_with_balanced{end};
        [r,~]=find(ismember(CP,i));
    end
end

save(['Results/concordant/' name '_concordant.mat'])


%% calculate full coupling (stoichiometric coupling) for comparison 
% network: metabolic network structure with fileds
network.stoichiometricMatrix = full(model.S);
network.reversibilityVector = (model.lb < 0 & model.ub > 0);
network.Reactions =  model.rxns;
network.Metabolites = model.mets;
 
coupling_matrix = F2C2('glpk',network);
 
save(['Results/stoichiometric_coupling/' name '_stoichiometric_coupling.mat'],'network','coupling_matrix')

%% get kinetic modules and generate final tables
%files_c = dir('Results/concordant/*.mat');
%fc = find(arrayfun(@(x) strcmp(files_c(x).name,[name '_concordant.mat']),1:size(files_c,1))); % index of file in folder Results/concordant 
fc = [name '_concordant.mat'];

cd Code/kinetic_modules/
disp('kinetic module')
system(['Rscript code_kineticModule_analysis.R' ' ' fc])


mkdir('../../Results/Overview/') 
mkdir('../../Results/MetSingle/') % absolute concentration robustness
mkdir('../../Results/MetDouble/') % absolute concentration ratio robustness
mkdir('../../Results/Reactions_Giant/') 
disp('Result tables')

fc = [name '_concordant.RData'];
system(['Rscript Get_excel_tables_kinetic_modules.R' ' ' fc])

end
