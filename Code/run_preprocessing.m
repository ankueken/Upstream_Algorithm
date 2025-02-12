function [model,name] = run_preprocessing(files,f,type,pathToR)

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

%% split into elementary reactions
disp('split reactions into elementary reactions')
if strcmp(type,'fixed') || strcmp(type,'ordered')
    model_elementary_fixed = split_into_elementary_rxns_v1(model,'fixed');
elseif strcmp(type,'random')
    model_elementary_random = split_into_elementary_rxns_v1(model,'random');
else 
    warning('model reactions are split assuming ordered binding, to assume random binding change input arguments')
end

if strcmp(type,'random')
    [sol.x,sol.f,sol.stat,sol.output]=linprog(-model_elementary_random.c,model_elementary_random.S(model_elementary_random.csense=='L',:),model_elementary_random.b(model_elementary_random.csense=='L'),model_elementary_random.S(model_elementary_random.csense=='E',:),model_elementary_random.b(model_elementary_random.csense=='E'),model_elementary_random.lb,model_elementary_random.ub);
    
    model_elementary_random = convertToIrreversible(model_elementary_random);
    % save file 
    mkdir('Models/temp/')
    save(strcat('Models/temp/',name,'_random.mat'),"model_elementary_random")
    
    cd Code/preprocessing/

    % pathToR = 'C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe';
    system(strjoin({pathToR, ' get_AY_matrix.r',strcat('../../Models/temp/',name,'_random.mat')}));
    
    cd ../../Models/temp/
    % random model
    load(strcat(name,'_random_A.dat'))
    load(strcat(name, '_random_complexes.mat'))
    model_elementary_random.A=eval(['spconvert(' strcat(name,'_random_A') ')']);
    model_elementary_random.complexes=complexes;
    load(strcat(name,'_random_Y.dat'))
    model_elementary_random.Y=eval(['spconvert(' strcat(name,'_random_Y') ')']);
    clear complexes *_random_A *_random_Y
    cd ../../
    
    save(['Models/models_with_elementary_steps/' name '_pre_balanced_random.mat'],'-v7.3')
else
    [sol.x,sol.f,sol.stat,sol.output]=linprog(-model_elementary_fixed.c,model_elementary_fixed.S(model_elementary_fixed.csense=='L',:),model_elementary_fixed.b(model_elementary_fixed.csense=='L'),model_elementary_fixed.S(model_elementary_fixed.csense=='E',:),model_elementary_fixed.b(model_elementary_fixed.csense=='E'),model_elementary_fixed.lb,model_elementary_fixed.ub);
    
    model_elementary_fixed = convertToIrreversible(model_elementary_fixed);
    mkdir('Models/temp/')
    save(strcat('Models/temp/',name,'_fixed.mat'),"model_elementary_fixed")
    
    cd Code/preprocessing/
    % pathToR = '"C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe"';
    system(strjoin({pathToR, ' get_AY_matrix.r',strcat('../../Models/temp/',name,'_fixed.mat')}));
    
    cd ../../Models/temp/
    % fixed model
    load(strcat(name,'_fixed_A.dat'))
    load(strcat(name, '_fixed_complexes.mat'))
    model_elementary_fixed.A=eval(['spconvert(' strcat(name,'_fixed_A') ')']);
    model_elementary_fixed.complexes=complexes;
    load(strcat(name,'_fixed_Y.dat'))
    model_elementary_fixed.Y=eval(['spconvert(' strcat(name,'_fixed_Y') ')']);
    clear complexes *_fixed_A *fixed_Y
    cd ../../
    
    save(['Models/models_with_elementary_steps/' name '_pre_balanced_fixed.mat'],'-v7.3')
end
system('rm -r Models/temp/')
if exist('model_elementary_fixed')
    model = model_elementary_fixed;
else
    model = model_elementary_random;
end

end
