function [At,CC] = pre_concordant_sampling(model,name,B,num_samples)

% Input: model  -   model struct
%        name   -   string specifying model name
%        B      -   index of balanced complexes
%        num_samples - number of samples per round
%
% Output: At    -   candidate pairs to check for concordance
%         CC    -   complex x complex matrix indicating relation
%                   balanced complexes (diagonal entry 1)
%                   concordant complexes (trivial -1, otherwise 2)

    CC = sparse(zeros(length(model.complexes)));
    CC(B,B) = 1;
    
    %% species degree
    species_degree_complexes = sum(model.Y~=0,2);
    % species_degree_reactions = sum(model.S~=0,2);
    
    % trivially concordant pairs
    idx = find(species_degree_complexes==2);
    for i=1:length(idx)
        tcc = find(model.Y(idx(i),:)~=0);
        if CC(tcc,tcc)~=1
            CC(tcc,tcc) = [0 -1;-1 0];
           % CC(logical(eye(size(CC)))) = 0;
        end
    end
    
    %% all possible pairs
    At=nchoosek(1:length(model.complexes),2);
    
    CC_temp=CC;
    CC_temp(B,:)=5;
    CC_temp(:,B)=5;
    
    [checked_rows,checked_cols]=find(triu(CC_temp,1)~=0);
    At = setdiff(At,[checked_rows checked_cols; checked_cols checked_rows],'rows');
    clear checked_rows checked_cols
    counter=0;
    At_old_size = size(At,1)*2;
    %%
    while size(At,1)<At_old_size*0.99 || counter<50
        At_old_size = size(At,1);
        sample_f=[];
        
        disp(name)
    
        options = optimset('linprog');
        options.Display = 'off';
    
        counter = counter+1;
        disp('start...')
    
        At_rows_to_check = unique(At(:,1));
    
        % while ~isempty(At_rows_to_check)
    
            fprintf('Round %i, done %i, At size %i \n',counter, length(At_rows_to_check),size(At,1))
    
            sv= model.lb + (model.ub-model.lb).*rand(length(model.lb), num_samples);
    
            for s=1:num_samples
    
                c = sparse([zeros(length(model.lb),1); model.ub-model.lb; model.ub-model.lb]);
                
                N = sparse([model.S zeros(size(model.S)) zeros(size(model.S));
                     eye(size(model.S,2)) eye(size(model.S,2)) -eye(size(model.S,2))]);
    
                b = sparse([model.b;sv(:,s)]);
    
                [X,mini,ExitFlag]=linprog(-c,[],[],N,b,[model.lb; -ones(length(model.lb)*2,1)*1000],[model.ub; ones(length(model.lb)*2,1)*1000],options);
    
                if ExitFlag==1
                    sample_f(:,s) = sparse(X(1:length(model.lb)));
                end
            end
            clear sv N b c X
	    fprintf('sampling %i done\n',num_samples)
                   sample=round(model.A*sample_f,9);
            clear sample_f
    
             [m,n]=size(sample);
                sample_reshaped = reshape(sample',1,n,m); % every row in sample
                candidates = std((sample+eps)./(sample_reshaped+eps),0,2,'omitnan')./mean((sample+eps)./(sample_reshaped+eps),2,'omitnan');
            disp(size(candidates))   
	     CC_temp = sparse(reshape(candidates,m,m)>0.01);
	     disp(size(CC_temp))
           clear samples samples_reshaped
       
        [checked_rows,checked_cols]=find(triu(CC_temp,1)~=0);
        At = setdiff(At,[checked_rows checked_cols; checked_cols checked_rows],'rows');
	    clear checked_rows checked_cols
    end
end