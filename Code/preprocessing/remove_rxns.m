function [model_field] = remove_rxns(model_field,NRXNS,inx)

    [nr,nc] = size(model_field);
    
    if nr == NRXNS
        model_field(inx,:) = [];
    elseif nc == NRXNS
        model_field(:,inx) = [];
    end

end