function model_elementary = split_into_elementary_rxns_v1(model,type)
% Function that split reactions (with assigned GPR rule) into elementary reaction steps.
% GPR rules are considered, for isoenzymes two enzymes are introduced
% enzyme complex formation not accounted, considered as one new enzyme
%
% x(1) | x(2) -> reaction catalyzed by E1 and duplicate catalyzed by E2
%                introduced
% x(1) & x(2) -> reaction catalyzed by E3, being the complex of E1 and E2,
%                is introduced
%
% INPUT:
% model - model struct, including at least stoichiometric matrix and GPR rules used
%         to split reactions
% type  - the splitting can be done in two ways
%          'fixed' - fixed order binding, only one order of binding is
%                    considered the order itself is choosen randomly
%          'random' - all possible orders of binding considered, see
%                     convenience kinetics
%
%% split reactions

% if the rule is empty but the enzyme is known add artificial gene to
% correctly split afterwards
if isfield(model,'rxnECNumbers')
inx = intersect(find(cellfun(@isempty,model.rxnECNumbers)==0), find(cellfun(@isempty,model.rules)));
model.rules(inx) = model.rxnECNumbers(inx);
end

% create a list of all possible enzymes, enzyme complexes
enzyme_list = cellfun(@(x) strsplit(x,{'|'}),model.rules,"UniformOutput",false);
enzyme_list = [enzyme_list{:}]';
enzyme_list = cellfun(@(x) strrep(x,' ',''), enzyme_list, 'UniformOutput',false);
enzyme_list = unique(enzyme_list);
enzyme_list(cellfun(@isempty,enzyme_list)) = []; % remove empty field

NMET = length(model.mets);
NENZ = length(enzyme_list);

% prepare matrix
model_elementary = struct();
model_elementary.S = sparse(nan(size(model.S,1)+length(enzyme_list),0));
model_elementary.b = [model.b;zeros(length(enzyme_list),1)];
model_elementary.mets = model.mets;
for temp=1:length(enzyme_list)
    model_elementary.mets{NMET+temp} = strcat('E_',num2str(temp));
end
model_elementary.metNames = [model.metNames; enzyme_list];
model_elementary.c = zeros(0,0);
model_elementary.rxns = cell(0,0);
model_elementary.rxnNames = cell(0,0);
model_elementary.lb = zeros(0,0);
model_elementary.ub = zeros(0,0);
model_elementary.unexp_rxn_inx = zeros(0,0);
model_elementary.rules = cell(0,0);
complex_list = zeros(size(model_elementary.S,1),0);

for i = 1:length(model.rxns)
% disp(i/length(model.rxns))
substrates = find(model.S(:,i)<0);
products = find(model.S(:,i)>0);
    % if reaction has no GPR rule
    if isempty(model.rules(i)) || isempty(model.rules{i}) || length(substrates)>4 || length(products)>4
        model_elementary.S(:,end+1) = [model.S(:,i); zeros(size(model_elementary.S,1)-size(model.S,1),1)];
        model_elementary.c(end+1,1) = model.c(i);
        model_elementary.lb(end+1,1) = model.lb(i);
        model_elementary.ub(end+1,1) = model.ub(i);
        model_elementary.rxns(end+1,1) = model.rxns(i);
        model_elementary.rxnNames(end+1,1) = model.rxnNames(i);
        model_elementary.unexp_rxn_inx(end+1,1) = i;
        model_elementary.rules(end+1,1) = model.rules(i);
    else
        if strcmp(type,'fixed')
            %%{
            % if reaction has GPR rule assigned
            substrates = find(model.S(:,i)<0);
            products = find(model.S(:,i)>0);
            [~,~,enzymes] = intersect(cellfun(@(x) strrep(x,' ',''),...
                strsplit(model.rules{i},'|'),'UniformOutput',false),enzyme_list); % remember inx in enzyme list

            for e_inx = 1:length(enzymes)
                %% formation of enzyme-substrate complex
                for s_inx = 1:length(substrates)

                    newcol = zeros(size(model_elementary.S,1),1);

                    complex_formed = zeros(size(complex_list,1),1);
                    complex_formed(substrates(1:s_inx)) = 1;
                    complex_formed(length(model.mets)+enzymes(e_inx)) = 1;

                    if s_inx==1 % substrate + enzyme
                        newcol([substrates(s_inx), length(model.mets)+enzymes(e_inx)]) = [model.S(substrates(s_inx),i), -1];
                    else % substrate plus previous complex
                        newcol([substrates(s_inx), find(model_elementary.S(:,end)>0)]) = [model.S(substrates(s_inx),i), -1];
                    end

                    model_elementary.S(:,end+1) = newcol; % substrates

                    % check if the complex is already in the model
                    if ~any(all(complex_list == complex_formed))
                        complex_list(:,end+1) = complex_formed;

                        model_elementary.S(end+1,end) = 1; % complex formed
                        model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(substrates(1:s_inx))),'_'), '_complex'];
                        model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(substrates(1:s_inx))),'_'), '_complex'];
                        model_elementary.b(end+1,1) = 0;
                    else
                        % find complex inx
                        c_inx = find(all(complex_list == complex_formed));
                        model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                    end

                    % set bounds
                    model_elementary.c(end+1,1) = model.c(i);
                    model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                    model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                    model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s' num2str(s_inx)];
                    model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz' num2str(enzymes(e_inx)) '_subst' num2str(s_inx)];
                    model_elementary.unexp_rxn_inx(end+1,1) = i;
                    model_elementary.rules(end+1,1) = model.rules(i);
                end % done with formation of substrate-enzyme complex

                %% formation of product-enzyme complex
                newcol = zeros(size(model_elementary.S,1),1);

                % substrate-enzyme complex from previous reaction
                if ~isempty(substrates)
                    newcol(find(model_elementary.S(:,end)>0)) = -1;
                else % we split exchange rxns like A <=> as A + E <=> AE <=> E
                    newcol(NMET+enzymes(e_inx),end) = -1;
                end
                model_elementary.S(:,end+1) = newcol;
                
                % product-enzyme complex we will obtain
                complex_formed = zeros(size(complex_list,1),1);
                complex_formed(products) = 1;
                complex_formed(length(model.mets)+enzymes(e_inx)) = 1;
                
                if ~isempty(products)
                % check if the complex is already in the model
                if ~any(all(complex_list == complex_formed)) 
                    complex_list(:,end+1) = complex_formed;

                    model_elementary.S(end+1,end) = 1; % complex formed
                    model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products)),'_'), '_complex'];
                    model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products)),'_'), '_complex'];
                    model_elementary.b(end+1,1) = 0;
                else
                    % find complex inx
                    c_inx = find(all(complex_list == complex_formed));
                    model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                end
                else
                    model_elementary.S(NMET+enzymes(e_inx),end) = 1;
                end

                model_elementary.c(end+1,1) = model.c(i,1);

                if model.lb(i) >= 0
                    model_elementary.lb(end+1,1) = 0;
                else
                    model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                end

                model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s_p_transition'];
                model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz_' num2str(enzymes(e_inx)) '_product_complex_formation'];
                model_elementary.unexp_rxn_inx(end+1,1) = i;
                model_elementary.rules(end+1,1) = model.rules(i);

                %% release of products from enzyme-product complex
                for p_inx = 1:length(products)
                    newcol = zeros(size(model_elementary.S,1),1);

                    if p_inx < length(products)
                        complex_formed = zeros(size(complex_list,1),1);
                        complex_formed(products((p_inx+1):length(products))) = 1;
                        complex_formed(length(model.mets)+enzymes(e_inx)) = 1;
                    end

                    % products positive value, complex negative
                    ec = find(model_elementary.S(:,end)>0);
                    ec = ec(find(contains(model_elementary.mets(ec),'complex')));
                    newcol([products(p_inx), ec]) = [model.S(products(p_inx),i), -1];

                    model_elementary.S(:,end+1) = newcol; % substrates

                    if p_inx<length(products)
                        % check if the complex is already in the model
                        if ~any(all(complex_list == complex_formed))
                            complex_list(:,end+1) = complex_formed;

                            model_elementary.S(end+1,end) = 1; % complex formed
                            model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products((p_inx+1):length(products)))),'_'), '_complex'];
                            model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products((p_inx+1):length(products)))),'_'), '_complex'];
                            model_elementary.b(end+1,1) = 0;
                        else
                            % find complex inx
                            c_inx = find(all(complex_list == complex_formed));
                            model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                        end
                    else
                        model_elementary.S(NMET+enzymes(e_inx),end) = 1;
                    end

                    % set bounds
                    model_elementary.c(end+1,1) = model.c(i);
                    model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                    model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                    model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_p' num2str(p_inx)];
                    model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz_' num2str(enzymes(e_inx)) '_product_' num2str(p_inx)];
                    model_elementary.unexp_rxn_inx(end+1,1) = i;
                    model_elementary.rules(end+1,1) = model.rules(i);
                end
            end % end loop isoenzymes
            %}
        else
            % if reaction has GPR rule assigned
            substrates = find(model.S(:,i)<0);
            products = find(model.S(:,i)>0);
            [~,~,enzymes] = intersect(cellfun(@(x) strrep(x,' ',''),...
                strsplit(model.rules{i},'|'),'UniformOutput',false),enzyme_list); % remember inx in enzyme list

            for e_inx = 1:length(enzymes)
                %% formation of enzyme-substrate complex
                for level = 1:length(substrates)
                    % for s_inx = 1:length(substrates)
                    if level==1
                        ls = substrates(nchoosek(1:length(substrates),1));

                        newcol = zeros(size(model_elementary.S,1),size(ls,1));

                        complex_formed = zeros(size(complex_list,1),size(ls,1));
                        complex_formed(ls,:) = eye(length(ls));
                        complex_formed(length(model.mets)+enzymes(e_inx),:) = 1;

                        newcol([ls', length(model.mets)+enzymes(e_inx)],:) = [diag(model.S(ls,i)); repmat(-1,1,length(ls))];

                        % check if the complex is already in the model
                        in_list = find(ismember(complex_formed',complex_list','rows'));
                        out_list = find(~ismember(complex_formed',complex_list','rows'));

                        if ~isempty(out_list)
                            for f=1:length(out_list)
                                model_elementary.S(1:size(newcol,1),end+1) = newcol(:,out_list(f)); % substrates

                                complex_list(:,end+1) = complex_formed(:,out_list(f));

                                model_elementary.S(end+1,end) = 1; % complex formed
                                model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                model_elementary.b(end+1,1) = 0;
                            end
                        end
                        if ~isempty(in_list)
                            for f=1:length(in_list)
                                % find complex inx
                                model_elementary.S(1:size(newcol,1),end+1) = newcol(:,in_list(f));

                                c_inx = find(all(complex_list == complex_formed(:,in_list(f))));
                                model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                            end
                        end

                        % set bounds
                        for f=1:size(ls,1)
                            model_elementary.c(end+1,1) = model.c(i);
                            model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                            model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                            model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s' num2str(ls(f,:))];
                            model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz' num2str(enzymes(e_inx)) '_subst' num2str(ls(f,:))];
                            model_elementary.unexp_rxn_inx(end+1,1) = i;
                            model_elementary.rules(end+1,1) = model.rules(i);
                        end
                        complex_formed_old = complex_formed;
                        % end % done with formation of substrate-enzyme complex
                    else
                        complex_formed_old_v2=zeros(size(complex_formed_old,1),0);
                        for c_inx=1:size(complex_formed_old,2)

                            ls = setdiff(substrates,find(complex_formed_old(:,c_inx)));
                            isec = intersect(substrates,find(complex_formed_old(:,c_inx)));

                            newcol = zeros(size(model_elementary.S,1),size(ls,1));

                            complex_formed = repmat(complex_formed_old(:,c_inx),1,size(ls,1));
                            complex_formed(ls,:) = eye(length(ls));
                            % complex_formed(length(model.mets)+enzymes(e_inx),:) = 1;

                            old_c_inx = find(strcmp(model_elementary.mets,['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(isec)),'_'), '_complex']));
                            newcol([ls', old_c_inx'],:) = [diag(model.S(ls,i)); repmat(-1,1,length(ls))];

                            % check if the complex is already in the model
                            in_list = find(ismember(complex_formed',complex_list','rows'));
                            out_list = find(~ismember(complex_formed',complex_list','rows'));

                            if ~isempty(out_list)
                                for f=1:length(out_list)
                                    model_elementary.S(1:size(newcol,1),end+1) = newcol(:,out_list(f)); % substrates

                                    complex_list(:,end+1) = complex_formed(:,out_list(f));

                                    model_elementary.S(end+1,end) = 1; % complex formed
                                    model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                    model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                    model_elementary.b(end+1,1) = 0;
                                end
                            end
                            if ~isempty(in_list)
                                for f=1:length(in_list)
                                    model_elementary.S(1:size(newcol,1),end+1) = newcol(:,in_list(f));

                                    % find complex inx
                                    co_inx = find(all(complex_list == complex_formed(:,in_list(f))));
                                    model_elementary.S(NMET+NENZ+co_inx,end) = 1;
                                end
                            end

                            % set bounds
                            for f=1:size(ls,1)
                                model_elementary.c(end+1,1) = model.c(i);
                                model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                                model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                                model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s' num2str(ls(f,:))];
                                model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz' num2str(enzymes(e_inx)) '_subst' num2str(ls(f,:))];
                                model_elementary.unexp_rxn_inx(end+1,1) = i;
                                model_elementary.rules(end+1,1) = model.rules(i);
                            end
                            complex_formed_old_v2 = [complex_formed_old_v2 complex_formed];
                        end
                        complex_formed_old = unique(complex_formed_old_v2','rows')';
                    end
                end

                %% formation of product-enzyme complex
                newcol = zeros(size(model_elementary.S,1),1);

                % substrate-enzyme complex from previous reaction
                 if ~isempty(substrates)
                    newcol(find(model_elementary.S(:,end)>0)) = -1;
                else % we split exchange rxns like A <=> as A + E <=> AE <=> E
                    newcol(NMET+enzymes(e_inx),end) = -1;
                end
                model_elementary.S(:,end+1) = newcol;

                

                % product-enzyme complex we will obtain
                complex_formed = zeros(size(complex_list,1),1);
                complex_formed(products) = 1;
                complex_formed(length(model.mets)+enzymes(e_inx)) = 1;

                if ~isempty(products)
                % check if the complex is already in the model
                if ~any(all(complex_list == complex_formed))
                    complex_list(:,end+1) = complex_formed;

                    model_elementary.S(end+1,end) = 1; % complex formed
                    model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products)),'_'), '_complex'];
                    model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(products)),'_'), '_complex'];
                    model_elementary.b(end+1,1) = 0;
                else
                    % find complex inx
                    c_inx = find(all(complex_list == complex_formed));
                    model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                end
                else
                    model_elementary.S(NMET+enzymes(e_inx),end) = 1;
                end

                model_elementary.c(end+1,1) = model.c(i,1);

                if model.lb(i) >= 0
                    model_elementary.lb(end+1,1) = 0;
                else
                    model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                end

                model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s_p_transition'];
                model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz_' num2str(enzymes(e_inx)) '_product_complex_formation'];
                model_elementary.unexp_rxn_inx(end+1,1) = i;
                model_elementary.rules(end+1,1) = model.rules(i);

                complex_formed_old = complex_formed;
                %% release of products from enzyme-product complex
                for level = 1:length(products)
                    if level==1
                        ls = products(nchoosek(1:length(products),1));

                        newcol = zeros(size(model_elementary.S,1),size(ls,1));

                        complex_formed = repmat(complex_formed_old,1,size(ls,1));
                        temp =  complex_formed(ls,:);
                        temp(logical(eye(size(ls,1))))=0;
                        complex_formed(ls,:) = temp; % remove one metabolite from complex

                        newcol(ls,:) = diag(model.S(ls,i)); % single product released
                        newcol(NMET+NENZ+find(ismember(complex_list',complex_formed_old','rows')),:) = -1; % complex used

                        in_list = find(ismember(complex_formed',complex_list','rows'));
                        out_list = find(~ismember(complex_formed',complex_list','rows'));
                        ec = find(all(complex_formed(1:NMET,:)==0));
                        in_list = setdiff(in_list,ec);
                        out_list = setdiff(out_list,ec);

                        if ~isempty(ec)
                            for f=1:length(ec)
                                model_elementary.S(1:size(newcol,1),end+1) = newcol(:,ec(f));

                                % find enzyme
                                model_elementary.S(find(complex_formed(:,ec(f))~=0),end) = 1;
                            end
                        end
                        if ~isempty(out_list)
                            for f=1:length(out_list)
                                model_elementary.S(1:size(newcol,1),end+1) = newcol(:,out_list(f)); % product released

                                complex_list(:,end+1) = complex_formed(:,out_list(f));

                                model_elementary.S(end+1,end) = 1; % complex formed
                                model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                model_elementary.b(end+1,1) = 0;
                            end
                        end
                        if ~isempty(in_list)
                            for f=1:length(in_list)
                                % find complex inx
                                model_elementary.S(1:size(newcol,1),end+1) = newcol(:,in_list(f));

                                c_inx = find(all(complex_list == complex_formed(:,in_list(f))));
                                model_elementary.S(NMET+NENZ+c_inx,end) = 1;
                            end
                        end

                        % set bounds
                        for f=1:size(ls,1)
                            model_elementary.c(end+1,1) = model.c(i);
                            model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                            model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                            model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s' num2str(ls(f,:))];
                            model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz' num2str(enzymes(e_inx)) '_prod' num2str(ls(f,:))];
                            model_elementary.unexp_rxn_inx(end+1,1) = i;
                            model_elementary.rules(end+1,1) = model.rules(i);
                        end
                        complex_formed_old = complex_formed;
                        % end % done with formation of substrate-enzyme complex

                    else
                        complex_formed_old_v2=zeros(size(complex_formed_old,1),0);
                        for c_inx=1:size(complex_formed_old,2)

                            ls = setdiff(products,find(complex_formed_old(:,c_inx)));
                            isec = intersect(products,find(complex_formed_old(:,c_inx)));

                            newcol = zeros(size(model_elementary.S,1),size(isec,1));

                            complex_formed = repmat(complex_formed_old(:,c_inx),1,size(isec,1));
                            temp =  complex_formed(isec,:);
                            temp(logical(eye(size(isec,1))))=0;
                            complex_formed(isec,:) = temp; % remove one metabolite from complex

                            old_c_inx = find(strcmp(model_elementary.mets,['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(isec)),'_'), '_complex']));
                            if length(isec)>1
                                newcol([isec', old_c_inx'],:) = [diag(model.S(isec,i)); repmat(-1,1,length(isec))];
                            else
                                newcol([isec', old_c_inx'],:) = [model.S(isec,i); -1];
                            end

                            % check if the complex is already in the model
                            in_list = find(ismember(complex_formed',complex_list','rows'));
                            out_list = find(~ismember(complex_formed',complex_list','rows'));

                            ec = find(all(complex_formed(1:NMET,:)==0));
                            in_list = setdiff(in_list,ec);
                            out_list = setdiff(out_list,ec);

                            if ~isempty(ec)
                                for f=1:length(ec)
                                    model_elementary.S(1:size(newcol,1),end+1) = newcol(:,ec(f));

                                    % find enzyme
                                    model_elementary.S(find(complex_formed(:,ec(f))~=0),end) = 1;
                                end
                            end
                            if ~isempty(out_list)
                                for f=1:length(out_list)
                                    model_elementary.S(1:size(newcol,1),end+1) = newcol(:,out_list(f)); % substrates

                                    complex_list(:,end+1) = complex_formed(:,out_list(f));

                                    model_elementary.S(end+1,end) = 1; % complex formed
                                    model_elementary.mets{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                    model_elementary.metNames{end+1,1} = ['E', num2str(enzymes(e_inx)), '_', strjoin(sort(model_elementary.mets(find(complex_formed(1:NMET,out_list(f))))),'_'), '_complex'];
                                    model_elementary.b(end+1,1) = 0;
                                end
                            end
                            if ~isempty(in_list)
                                for f=1:length(in_list)
                                    model_elementary.S(1:size(newcol,1),end+1) = newcol(:,in_list(f));

                                    % find complex inx
                                    co_inx = find(all(complex_list == complex_formed(:,in_list(f))));
                                    model_elementary.S(NMET+NENZ+co_inx,end) = 1;
                                end
                            end

                            % set bounds
                            for f=1:size(ls,1)
                                model_elementary.c(end+1,1) = model.c(i);
                                model_elementary.lb(end+1,1) = min([-1000,min(model.lb)]);
                                model_elementary.ub(end+1,1) = max([1000,max(model.ub)]);
                                model_elementary.rxns{end+1,1} = [model.rxns{i} '_e' num2str(enzymes(e_inx)) '_s' num2str(ls(f,:))];
                                model_elementary.rxnNames{end+1,1} = [model.rxnNames{i} '_enz' num2str(enzymes(e_inx)) '_prod' num2str(ls(f,:))];
                                model_elementary.unexp_rxn_inx(end+1,1) = i;
                                model_elementary.rules(end+1,1) = model.rules(i);
                            end
                            complex_formed_old_v2 = [complex_formed_old_v2 complex_formed];
                        end
                        complex_formed_old = unique(complex_formed_old_v2','rows')';
                    end
                end 
            end 
        end% choose type
        model_elementary.csense = repmat('E',size(model_elementary.mets));
    end
end % end loop rxns
end