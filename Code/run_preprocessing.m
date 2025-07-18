function [model,name] = run_preprocessing(files,f)

% the code runs the preprocessing of models in folder Models/original 
% and split of reactions into elementary reaction steps under 
% the assumption of a fixed ordered binding
%
% Input: files - result from dir() where model files are located (e.g. files = dir('Model/original/*.xml');)
%        f     - (integer) index of the model file for which the code should be executed
%        type  - (string) the splitting can be done in two ways
%                'fixed' - fixed order binding, only one order of binding is
%                    considered the order itself is choosen randomly
%                'random' - all possible orders of binding considered, see
%                    convenience kinetics
%       pathToR - absolute path to R (e.g. pathToR = '"C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe"';))
% Output: updated model file also written into Models/models_with_elementary_steps/
%
% Requirements: CobraToolbox, R with packages igraph and R.matlab
% R will be called by systems('pathToR R_file.r ags')


disp(f)
addpath(genpath('Code/'))
%%{
if endsWith(files(f).name,'.sbml') 
    name = files(f).name(1:end-5)
else
    name = files(f).name(1:end-4)
end

% disp('set cobra path')
% COBRA_PATH = '/work/ankueken/Git/cobratoolbox/';
% addpath(genpath(COBRA_PATH));

model = readCbModel(strcat(files(f).folder,'/',files(f).name));

disp('Done read model')

%% clean model
% we want to use generic bounds and only account for reversibility
model.lb(model.lb<0) = -1000;
model.lb(model.lb>0) = 0;
model.ub(model.ub<0) = 0;
model.ub(model.ub>0) = 1000;

model.lb(model.c~=0) = 0;
model.ub(model.c~=0) = 1000;

if ~isfield(model,'csense')
    model.csense = repmat('E',size(model.mets));
end

disp('Clean model')
model=removeRxns(model,model.rxns(find(all(model.S==0))));
model=removeMetabolites(model,model.mets(find(all(model.S'==0))));

[solo.x,solo.f,solo.stat,solo.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub);

[mini,maxi] = linprog_FVA(model,0.001);
thr=1e-9;
BLK=model.rxns(find(abs(mini)<thr & abs(maxi)<thr));

model=removeRxns(model,BLK);

[sol.x,sol.f,sol.stat,sol.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub);

if abs(sol.f)<abs(solo.f)*0.5
    disp('No biomass due to removal of blocked reactions')
    % save(['Results/Problems/' name '.mat'])
    return
end

clear BLK sol solo files

model_elementary = convertToIrreversible(model);
mkdir('Models/temp/')
save(strcat('Models/temp/',name,'.mat'),"model_elementary")
    
    cd Code/preprocessing/
    % pathToR = '"C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe"';
    system(strjoin({'Rscript get_AY_matrix.r',strcat('../../Models/temp/',name,'.mat')}));
    
    cd ../../Models/temp/
    % fixed model
    load(strcat(name,'_A.dat'))
    load(strcat(name, '_complexes.mat'))
    model_elementary.A=eval(['spconvert(' strcat(name,'_A') ')']);
    model_elementary.complexes=complexes;
    load(strcat(name,'_Y.dat'))
    model_elementary.Y=eval(['spconvert(' strcat(name,'_Y') ')']);
    clear complexes *_A *_Y
    cd ../../
    
    save(['Models/models_irrev/' name '_pre_balanced.mat'],'-v7.3')

system('rm -r Models/temp/')
model = model_elementary;

end
