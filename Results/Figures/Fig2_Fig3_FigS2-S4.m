T=readtable("Table_ED-1.xlsx",NumHeaderLines=3,ReadVariableNames=true);

%% Figure 2
bar(log2([T.maxSize(strcmp(T.BindingType,'ordered')),T.maxSize(strcmp(T.BindingType,'random'))]),'grouped')
ylabel('Number of complexes in giant module')
legend('ordered binding','random binding')
legend boxoff
labels={
    '{\itA. baumannii} iCN718'
    '{\itA. thaliana} AraCore'
    '{\itA. thaliana} iRS1597'
    '{\itB. subtilis} iYO844'
    '{\itC. variabilis} iAJ526'
    '{\itC. difficile} iCN900'
    '{\itC. ljungdahlii} iHN637'
    '{\itE. coli} iAF1260b'
    '{\itE. coli} iJR904'
    '{\itH. pylori} iIT341'
    '{\itH. sapiens} iAB_RBC_283'
    '{\itL. lactis} iNF517'
    '{\itM. barkeri} iAF692'
    '{\itM. tuberculosis} iNJ661'
    '{\itP. berghei} iAM_Pb448'
    '{\itP. cynomolgi} iAM_Pc455'
    '{\itP. falciparum} iAM_Pf480'
    '{\itP. knowlesi} iAM_Pk459'
    '{\itP. vivax} iAM_Pv461'
    '{\itP. putida} iJN746'
    '{\itS. cerevisiae} iMM904'
    '{\itS. cerevisiae} iND750'
    '{\itS. dysenteriae} iSDY_1059'
    '{\itS. flexneri} iS_1188'
    '{\itS. aureus} iSB619'
    '{\itS. aureus} iYS854'
    '{\itS. elongatus} iJB785'
    '{\itSynechocystis sp.} iJN678'
    '{\itSynechocystis sp.} iSynCJ816'
    '{\itT. maritima} iLJ478'
    '{\itT. cruzi} iIS312 Amastigote'
    '{\itT. cruzi} iIS312 Epimastigote'
    '{\itT. cruzi} iIS312 Trypomastigote'
    '{\itT. cruzi} iIS312'
};
labels=strrep(labels,'_','-');
set(gca,'XTick',1:height(T)/2,'XTickLabel',labels,'XTickLabelRotation',45,'YTickLabel',2.^yticks)

%% Figure ED-2

% ordered giant - reactions
subplot(2,2,1)
plot(T.x_reactions(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')),'.k')
[x,p]=corr(T.x_reactions(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')))
ylabel('Size of giant module')
text(T.x_reactions(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')),num2str([1:34]'),'Color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')
% ordered giant - metabolites
subplot(2,2,2)
plot(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')),'.k')
[x,p]=corr(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')))
text(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'ordered')), T.maxSize(strcmp(T.BindingType,'ordered')),num2str([1:34]'),'Color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')
% random giant - reactions
subplot(2,2,3)
plot(T.x_reactions(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')),'.k')
[x,p]=corr(T.x_reactions(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')))
ylabel('Size of giant module')
xlabel('Number of reactions')
text(T.x_reactions(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')),num2str([1:34]'),'Color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')% random giant - metabolites
subplot(2,2,4)
plot(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')),'.k')
[x,p]=corr(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')))
xlabel('Number of metabolites')
text(T.x_metabolites_EnzymesAndEnzyme_metaboliteComplexes(strcmp(T.BindingType,'random')), T.maxSize(strcmp(T.BindingType,'random')),num2str([1:34]'),'Color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','center')% random giant - metabolites

%% Figure ED-3
mean(T.x_ReactionsWithSubstrateComplexInGiantComponent_1(find(strcmp(T.BindingType,'ordered') & T.maxSize>9)))

subplot(2,1,1)
bar(T.x_ReactionsWithSubstrateComplexInGiantComponent_1(strcmp(T.BindingType,'ordered')))
ylabel('Percentage of reactions with substrate complex in giant module')
subplot(2,1,2)
bar(T.x_ReactionsWithSubstrateComplexInGiantComponent_1(strcmp(T.BindingType,'random')))
ylabel('Percentage of reactions with substrate complex in giant module')
set(gca,'XTick',1:height(T)/2,'XTickLabel',labels,'XTickLabelRotation',45)

%% Figure 3
figure
subplot(1,2,1)
barh(flip(log2(T.x_DistinctFreeMetabolites_notBoundToEnzyme_(strcmp(T.BindingType,'ordered')))),'FaceColor',[.5 .5 .5])
xlim([0 8])
set(gca,'XTick',1:2:max(xticks))
set(gca,'XTickLabels',2.^xticks,'XDir','reverse','YAxisLocation','right','YTickLabel',[],'FontSize',12)

xlabel({'Number of metabolites with';'absolute concentration robustness'})
txt=labels;
txt=flip(txt);
for i=1:length(txt)
    text(-3,i,txt(i), 'HorizontalAlignment', 'center','FontSize',10)
end
subplot(1,2,2)
barh(flip(log2(T.x_DistinctFreeMetabolitePairs_metabolitesNotBoundToEnzyme_(strcmp(T.BindingType,'ordered')))),'FaceColor',[.5 .5 .5])
xlim([0 8])
set(gca,'XTick',1:2:max(xticks))
set(gca,'XTickLabels',2.^xticks,'YTickLabel',[],'OuterPosition',[0.56 0 0.4097 1],'FontSize',12)
xlabel({'Number of metabolite pairs with';'absolute concentration ratio robustness'})

%% Figure ED-4

bar(log2([T.maxSize(strcmp(T.BindingType,'ordered')),T.maxSize_1(strcmp(T.BindingType,'ordered'))]),'grouped')
ylabel('Number of complexes in giant module')
legend('kinetic coupling','stoichiometric coupling')
legend boxoff
set(gca,'XTick',1:height(T)/2,'XTickLabel',labels,'XTickLabelRotation',45,'YTickLabel',2.^yticks)


