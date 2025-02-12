clear
OrganismTable = readtable("../../../../Book1.xlsx","NumHeaderLines",5);
files=dir('Overview*.csv');

BindingType = cell(length(files),1);
Model = cell(length(files),1);
Organism = cell(length(files),1);
T = zeros(length(files),26);

for f=1:length(files)
    Trow = readtable(files(f).name);

    %% general info
    if contains(files(f).name,'random')
        BindingType{f,1} = 'random';
    else
        BindingType{f,1} = 'ordered';
    end

    Model(f,1) = {strrep(strrep(strrep(files(f).name,'Overview_',''),...
        '_concordant_random.csv',''),...
        '_concordant_fixed.csv','')};

    % Organism(f,1) = OrganismTable.Var2(find(contains(OrganismTable.Var1,[Model{f,1} '.xml'])));

    % complexes, reactions, free metabolites, metabolites, enzymes and enzyme-metabolite complexes
    T(f,[1 2 4]) = Trow.Var2([1 2 3]);

    % if contains(files(f).name,'fixed')
        load(['../../concordant_final/' strrep(strrep(files(f).name,'Overview_',''),'.csv','.mat')],'Results_balanced')
    % else
    %     load(['../concordant_random/' strrep(strrep(files(f).name,'Overview_',''),'.csv','.mat')],'Results_balanced')
    % end
    enzyme_inx = find(strcmp(Results_balanced.MODEL_r{1}.mets,'E_1'));
    T(f,3) = enzyme_inx-1; % number of free metabolites

    %% kinetic coupling
    T(f,5:9) = Trow.Var2(4:8);
    T(f,10) = (T(f,7)/T(f,1))*100;
    T(f,11) = Trow.Var2(23);
    T(f,12) = (T(f,11)/T(f,2))*100;

    %% robustness of concentrations (ratios)
    % singletons + pairs
    T(f,[13 15]) = Trow.Var2([14 15]);
    T(f,14) = (T(f,13)/T(f,4))*100;
    T(f,16) = (T(f,15)/(T(f,4)*(T(f,4)-1)/2))*100;

    MetSingle=readtable(['../MetSingle/MetSingle_' strrep(files(f).name,'Overview_','')],'ReadVariableNames',true,'Delimiter','\t');
    MetSingle=table2cell(MetSingle);
    MetSingle=cellfun(@(x) strrep(x,' "',''),MetSingle,'UniformOutput',false);
    MetSingle=cellfun(@(x) strrep(x,'"',''),MetSingle,'UniformOutput',false);
    
    M_MetSingle = intersect(MetSingle,...
    Results_balanced.MODEL_r{1}.mets(1:enzyme_inx-1));
    
    M_MetSingle_length = length(M_MetSingle);
    writecell(M_MetSingle,['../MetSingle/M_MetSingle_' strrep(files(f).name,'Overview_','')])
    T(f,17) = M_MetSingle_length;
    T(f,18) = (T(f,18)/T(f,3))*100;

    MetDouble=readtable(['../MetDouble/MetDouble_' strrep(files(f).name,'Overview_','')],'ReadVariableNames',false,'Delimiter','"');
    if size(MetDouble)==[2 2]
	MetDouble=rows2vars(MetDouble);
	MetDouble=MetDouble(2,2:end);
    end

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
   
    [~,M_MetDouble1] = intersect(MetDouble(:,1),...
    Results_balanced.MODEL_r{1}.mets(1:enzyme_inx-1));

    [~,M_MetDouble2] = intersect(MetDouble(:,2),...
    Results_balanced.MODEL_r{1}.mets(1:enzyme_inx-1));

    M_MetDouble = MetDouble(intersect(M_MetDouble1,M_MetDouble2),:);

    M_MetDouble_length = size(M_MetDouble,1);
    writecell(M_MetDouble,['../MetDouble/M_MetDouble_' strrep(files(f).name,'Overview_','')])
    T(f,19) = M_MetDouble_length;
    T(f,20) = (T(f,19)/(T(f,3)*(T(f,3)-1)/2))*100;

    %% clustering of coupled metabolites
    T(f,21:25) = Trow.Var2(16:20);
    T(f,26) = (T(f,23)/T(f,4))*100;

end

Results = table(BindingType,Model,Organism,T);
writetable(Results,'Results_table.xlsx')