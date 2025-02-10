function [model,BLK] = clean_model(name,model,constraints,pathToR,threshold)
% function to preprocess models used to identify balanced complexes
%
% [model,BLK] = clean_model(name,model,constraints,threshold,pathToR,BLK)
%
% Input
%       name: name used to save model files
%       model: struct with at least following fields
%           .c      objective vectore
%           .S      stoichiometric matrix
%           .b      right-hand side vector
%           .lb     lower bound on flux
%           .ub     upper bound on flux
%           .rxns   cell array of reaction names
%           .mets   cell array of metabolite names
%       constraints: specify scenario
%                    'all_reversible' - all reactions considered reversible
%                    'irreversibility_considered' - reaction irreversibility considered
%                    'objective' - maximize objective function 
%                                  (i.e. biomass)
%       pathToR: specify path to R in the form given below 
%                e.g. '"C:\Program Files\R\R-3.6.1\bin\Rscript.exe"'
%       threshold: value when to consider complex balanced (e.g. 1e-9)
%                  (optional)
%
% Output
%       model: preprocessed model, extended by fields .A and .Y
%              no blocked reactions, reversible reactions split
%              constraints according to specified scenario
%              ('all_reversible', 'irreversibility_considered', 'objective')
%       BLK: indices of reactions blocked in the original model
%
% Step 1:
%   check if there are zero rows or columns in the stoichiometric matrix
%   and remove them
% Step 2:
%   set case (irreversibility considered, objective)
%   and convert to irreversibility
% Step 3:
%   identify and remove blocked reactions
% Step 4:
%   calculate A and Y matrix
%
% -----------------------------------------------------------
if nargin<6
    threshold=1e-9;
end

% Step 1. Remove zero rows / columns
model=removeRxns(model,model.rxns(find(all(model.S==0))));
model=removeMetabolites(model,model.mets(find(all(model.S'==0))));

% Step 2. Set constraints

% find biomass coefficient
if all(model.c==0)
    obj=find(contains(model.rxns,'bio','IgnoreCase',true));
    % N. pharaonis: model.c(200)=1; find(strcmp(model.rxns,'RM00001'))
    if isempty(obj)
        disp('Please specify flux to maximize in model.c')
        keyboard
    else
        model.c(obj(1))=1;
    end
end

% convert to irreversible (split reversible reactions)
model.lb(model.c~=0)=0; % biomass reaction always considered irreversible
model=convertToIrreversible(model);

% Find and remove blocked reactions
[mini,maxi]=linprog_FVA(model,0);
BLK=model.rxns(intersect(find(abs(mini)<threshold),find(abs(maxi)<threshold)));
model=removeRxns(model,BLK);
[~,Sol.f]=linprog(-model.c,[],[],model.S,model.b,model.lb,model.ub);

% reaction irreversibility considered
if strcmp(constraints,'irreversibility_considered')
    model.c(:)=0;
    [~,Sol.f]=linprog(-model.c,[],[],model.S,model.b,model.lb,model.ub);
end

% Set objective to maximize
if strcmp(constraints,'objective')
    [~,Sol.f]=linprog(-model.c,[],[],model.S,model.b,model.lb,model.ub);
    model.lb(model.c~=0)=Sol.f*-0.99;
    model.ub(model.c~=0)=Sol.f*-0.991;
end

% Step 4. Compute A and Y matrix and save model

if strcmp(constraints,'all_reversible')
    save(strcat('Results/final_models_all_reversible/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_all_reversible/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_all_reversible/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_all_reversible/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_all_reversible/',name,'_final_model.mat'))
    
elseif strcmp(constraints,'irreversibility_considered')
    save(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_irreversibility_considered/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_irreversibility_considered/',name,'_final_model.mat'))
    
elseif strcmp(constraints,'objective')
    save(strcat('Results/final_models_objective/',name,'_final_model.mat'))
    cd functions
    system(strjoin({pathToR, 'get_AY_matrix.r',strcat('../Results/final_models_objective/',name,'_final_model.mat')}));
    cd ..
    A=importdata(strcat('Results/final_models_objective/',name,'_final_model.A'));
    model.A=A.data;
    model.complexes=A.textdata(2:end);
    Y=importdata(strcat('Results/final_models_objective/',name,'_final_model.Y'));
    model.Y=Y.data;
    clear A Y
    save(strcat('Results/final_models_objective/',name,'_final_model.mat'))
end

end