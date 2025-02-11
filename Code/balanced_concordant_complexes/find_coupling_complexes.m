function CC = find_coupling_complexes(model,group,CC,At,start,stop)

options = optimset('linprog');
options.Display = 'off';

%% setting start and stop allows to split the candidate set into several jobs
% if start=1 and stop=length(At(:,2)) we run the whole set of candidate pairs 
if ~exist('stop') || stop>length(At(:,2))
    stop=length(At(:,2));
end
if ~exist('start')
    start=1;
end
if start>length(At(:,2))
    disp('Done')
    return
end

%%
disp(length(At(:,2)))
disp('start...')

a = sparse([eye(size(model.S,2)) -model.ub; -eye(size(model.S,2)) model.lb]);

for i=start:stop
    % disp(i)

    if At(i,1)~=At(i,2) && group(At(i,2))~='N'

        % maximize

        [R.x,R.f_k,R.ExitFlag]=linprog([-model.A(At(i,1),:) 1],a,[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At(i,2),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; 0],[ones(size(model.S,2),1)*1e9; 999],options);

        if R.ExitFlag == 1

            Maximum_c_p = (model.A(At(i,1),:)*R.x(1:end-1))/(model.A(At(i,2),:)*R.x(1:end-1));

            % minimize
            [R.x,R.f_k,R.ExitFlag]=linprog([model.A(At(i,1),:) 1],a,[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At(i,2),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; 0],[ones(size(model.S,2),1)*1e9; 999],options);

            if R.ExitFlag == 1

                Minimum_c_p = (model.A(At(i,1),:)*R.x(1:end-1))/(model.A(At(i,2),:)*R.x(1:end-1));
            else
                Minimum_c_p = -Inf;
            end
        else
            Minimum_c_p = -Inf; Maximum_c_p = Inf;
        end
    end
    if At(i,1)~=At(i,2) && group(At(i,2))~='P'

        [R.x,R.f_k,R.ExitFlag]=linprog([-model.A(At(i,1),:) 1],a,[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At(i,2),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; -999],[ones(size(model.S,2),1)*1e9; 0],options);

        if R.ExitFlag == 1

            Maximum_c_n = (model.A(At(i,1),:)*R.x(1:end-1))/(model.A(At(i,2),:)*R.x(1:end-1));

            [R.x,R.f_k,R.ExitFlag]=linprog([model.A(At(i,1),:) 1],a,[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At(i,2),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; -999],[ones(size(model.S,2),1)*1e9; 0],options);

            if R.ExitFlag == 1

                Minimum_c_n = (model.A(At(i,1),:)*R.x(1:end-1))/(model.A(At(i,2),:)*R.x(1:end-1));
            else
                Minimum_c_n = -Inf;
            end

        else
            Maximum_c_n = Inf; Minimum_c_n = -Inf;
        end
    end
    % update CC, concordant pair denoted by value 2
    if group(At(i,2)) == 'P' && CC(At(i,1),At(i,2)) == 0
        CC(At(i,1),At(i,2)) = (round(Maximum_c_p,2) == round(Minimum_c_p,2))*2;
    elseif group(At(i,2)) == 'N' && CC(At(i,1),At(i,2)) == 0
        CC(At(i,1),At(i,2)) = (round(Maximum_c_n,2) == round(Minimum_c_n,2))*2;
    elseif CC(At(i,1),At(i,2)) == 0
        CC(At(i,1),At(i,2)) = (round(Maximum_c_p,2) == round(Minimum_c_p,2) & ...
            round(Maximum_c_n,2) == round(Minimum_c_n,2) & ...
            round(Maximum_c_p,2) == round(Minimum_c_n,2))*2;
    end
end

% save(['Results/concordant_random/' strrep(name,'_pre','') '_' num2str(size(At,1)) '_' num2str(start) '_' num2str(stop) '.mat'],'-v7.3')
% save(['Results/concordant_fixed/' strrep(name,'_pre','') '_' num2str(size(At,1)) '_' num2str(start) '_' num2str(stop) '.mat'],'-v7.3')

end
