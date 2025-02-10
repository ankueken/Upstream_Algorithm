function [B,group] = find_balanced_complexes(model,threshold)
% function to find balanced complexes
%
% B = find_balanced_complexes(model,threshold)
%
% Input
%   model: struct with at least following fields
%           .c      objective vectore
%           .S      stoichiometric matrix
%           .b      right-hand side vector
%           .A      complex-reaction matrix
%           .lb     lower bound on flux
%           .ub     upper bound on flux
%   threshold: value when to consider complex balanced (e.g. 1e-9)
%              (optional)
% 
% Output
%   B:  vector of indices of balanced complexes


if nargin<2
    threshold=1e-9;
end

options = optimset('linprog');
options.Display = 'off';

for i=1:size(model.A,1)

    model.c = model.A(i,:)';
    
    [~,Sol_max.f,Sol_max.stat]=linprog(model.c*-1,[],[],model.S,model.b,model.lb,model.ub,options);
    if Sol_max.stat == 1
        Maximum(i) = Sol_max.f*-1;
    else
        Maximum(i) = nan;
    end
    
    [~,Sol_min.f,Sol_min.stat]=linprog(model.c,[],[],model.S,model.b,model.lb,model.ub,options);
    if Sol_min.stat ==1
        Minimum(i) = Sol_min.f;
    else
        Minimum(i) = nan;
    end
end

Minimum=round(Minimum*(1/threshold))/(1/threshold);
Maximum=round(Maximum*(1/threshold))/(1/threshold);

B = intersect(find(Maximum==0),find(Minimum==0));

group = repmat('U',length(Maximum),1);
group(intersect(find(Maximum==0),find(Minimum==0))) = 'B';
group(intersect(find(Maximum>0),find(Minimum>=0))) = 'P';
group(intersect(find(Maximum<=0),find(Minimum<0))) = 'N';

end