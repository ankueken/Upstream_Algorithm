% the script takes out the enzymes and enzyme-substrate complexes found 
% the free metabolites are then written to M_MetSingle

files=dir('concordant/*.mat')
T=table([],[],[],'VariableNames',{'model','M_MetSingle','M_total'});
for f=1:length(files)

name = files(f).name(1:end-4);

try
    MetSingle=readtable(['MetSingle/MetSingle_' name '.csv'],'ReadVariableNames',true,'Delimiter','\t');
    MetSingle=table2cell(MetSingle);
    MetSingle=cellfun(@(x) strrep(x,' "',''),MetSingle,'UniformOutput',false);
    MetSingle=cellfun(@(x) strrep(x,'"',''),MetSingle,'UniformOutput',false);


    load([files(f).folder '/' files(f).name],'Results_balanced')

    enzyme_inx = find(strcmp(Results_balanced.MODEL_r{1}.mets,'E_1'));
    
    M_MetSingle = intersect(MetSingle,...
    Results_balanced.MODEL_r{1}.metNames(1:enzyme_inx-1));
    
    M_MetSingle_length = length(M_MetSingle);
    writecell(M_MetSingle,['MetSingle/M_MetSingle_' name '.csv'])
    T = [T; table({name},M_MetSingle_length,enzyme_inx-1,'VariableNames',{'model','M_MetSingle','M_total'})];
catch
    
end

end


files=dir('concordant/*.mat')
T=table([],[],[],'VariableNames',{'model','M_MetDouble','M_total'});
for f=1:length(files)

name = files(f).name(1:end-4);

try
    MetDouble=readtable(['MetDouble/MetDouble_' name '.csv'],'ReadVariableNames',false,'Delimiter','"');
    MetDouble=table2cell(MetDouble);
    keep=[];
    for c=1:size(MetDouble,2)
	try
	NN=all(cellfun(@isnan,MetDouble(:,c))==1);
	catch
	keep = [keep,c];
	end
    end
    MetDouble = MetDouble(:,keep);

 
    load([files(f).folder '/' files(f).name],'Results_balanced')

    enzyme_inx = find(strcmp(Results_balanced.MODEL_r{1}.mets,'E_1'));
    
    [~,M_MetDouble1] = intersect(MetDouble(:,1),...
    Results_balanced.MODEL_r{1}.metNames(1:enzyme_inx-1));

    [~,M_MetDouble2] = intersect(MetDouble(:,2),...
    Results_balanced.MODEL_r{1}.metNames(1:enzyme_inx-1));

    M_MetDouble = MetDouble(intersect(M_MetDouble1,M_MetDouble2),:);

    
    M_MetDouble_length = size(M_MetDouble,1);
    writecell(M_MetDouble,['MetDOuble/M_MetDouble_' name '.csv'])
    T = [T; table({name},M_MetDouble_length,enzyme_inx-1,'VariableNames',{'model','M_MetDouble','M_total'})];
catch
    
end

end




