function BLK = linprog_FindBlockedRxns(model,opt_percent,thr)
% opt_percent value between 0 and 1
if ~exist("opt_percent")
    opt_percent=1;
end
if ~exist("thr")
    thr=1e-9;
end

options = optimset('linprog');
options.Display = 'off';

model.lb(model.c~=0)=0;
model.ub(model.c~=0)=1000;

% find optimium biomass
[Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub,options);
% fix it to range
model.lb(model.c~=0)=Sol.f*-opt_percent;
model.ub(model.c~=0)=Sol.f*-1;

for i=1:size(model.S,2)
    model.c(:)=0;
    model.c(i)=1;
    [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub,options);
    if Sol.stat~=1
        maxi(i,1)=0;
    else
        maxi(i,1)=Sol.f*-1;
    end
end
to_test = find(maxi<thr);
for i=1:length(to_test)
    model.c(:)=0;
    model.c(to_test(i))=1;
    [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub,options);
    if Sol.stat~=1
        mini(to_test(i),1)=0;
    else
        mini(to_test(i),1)=Sol.f;
    end
end

BLK=model.rxns(intersect(find(abs(mini)<thr),find(abs(maxi)<thr)));
end