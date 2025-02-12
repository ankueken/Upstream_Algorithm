clear 
load concordant_final\kinetic_Ecoli_GEM_Maranas.mat
model=Results_balanced.MODEL_r{1};

model.c(4552) = 1;

sol_opt=optimizeCbModel(model);

for i=1:length(model.rxns)
    model_e = model;
    model_e.lb(i)=0;
    model_e.ub(i)=0;
    sol_e=optimizeCbModel(model_e);

    bio_KO(i,1) = sol_e.f;
end

Rxns_giant=readtable('Reactions_Giant\Reactions_Giant_kinetic_Ecoli_GEM_Maranas.csv')

essential_giant = sum(bio_KO(Rxns_giant.Var2)<1e-3);
non_essential_giant = sum(bio_KO(Rxns_giant.Var2)>1e-3);

essential_net = sum(bio_KO(setdiff(1:length(model.rxns),Rxns_giant.Var2))==0);
non_essential_net = sum(bio_KO(setdiff(1:length(model.rxns),Rxns_giant.Var2))>0);

FT_Table = [essential_giant essential_net; non_essential_giant non_essential_net];
[H,p]=fishertest(FT_Table);

FT_Table(1,:)./FT_Table(2,:)