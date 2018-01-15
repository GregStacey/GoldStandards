
% also see /Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Justification document/justifyNoRef.m
% this ^ includes analysis of our replicates


%% Read in allComplexes.csv

fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Input/allComplexes.csv';
fid = fopen(fn);
head = fgetl(fid);
head = strsplit(head,';');
head = head(1:12);

corum = cell(10000,size(head,2));
corumN = nan(10000,1);
cc = 0;
while ~feof(fid)
  cc = cc+1;
  t = fgetl(fid);
  t1 = strsplit(t,';','CollapseDelimiters',0);
  corum(cc,:) = t1;
  corumN(cc) = length(strsplit(t1{5},','));
end
corum = corum(1:cc,:);
corumN = corumN(1:cc);




%% Now do it for all datasets

datasets = {'Anders' 'Nick_HeLa' 'Nick_Ap' 'Nick_tiss' 'Jenny' 'Craig' 'Mike' 'Mike_{10}'};
maindir = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Output/tmp/';
fn = dir([maindir 'dataset*']);

fbad = nan(size(fn));
fbad_corum = nan(length(fn),size(corum,1));
ravg_corum = nan(length(fn),size(corum,1));
eluteTime_corum = nan(length(fn),size(corum,1));
eluteTimestd_corum = nan(length(fn),size(corum,1));
ds = nan(size(fn));
for ii = 1:length(fn)
  load([maindir fn(ii).name])
  
  ds(ii) = str2double(fn(ii).name(8));
  
  RR = -1*Dist.R;
  a = triu(TP_Matrix);
  I1 = a(:)==1;
  fbad(ii) = sum(RR(I1)<0) / sum(I1)
  
  [~,Ipeak] = max(Chromatograms,[],2);
  
  % find the average R, average elution time for each corum complex
  for jj = 1:size(corum,1)
    s1 = corum{jj,5};
    s1 = strrep(s1,')','');
    s1 = strrep(s1,'(','');
    prots = strsplit(s1,',');
    I = [];
    for kk = 1:length(prots)
      I = [I; find(~cellfun('isempty',strfind(Protein.NoIsoform,prots{kk})))];
    end
    if isempty(I);continue;end
    
    tmp = RR(I,I);
    tmp(tmp==1) = nan;
    N = length(tmp)^2 - length(tmp);
    fbad_corum(ii,jj) = sum(tmp(:)<0) / N;
    ravg_corum(ii,jj) = nanmean(tmp(:));
    eluteTime_corum(ii,jj) = nanmean(Ipeak(I));
    eluteTimestd_corum(ii,jj) = nanstd(Ipeak(I));
  end
end

DS = unique(ds);
avg_fbad = nan(size(DS));
for ii = 1:length(DS)
  avg_fbad(ii) = nanmean(fbad(ds==DS(ii) & fbad~=0));
end



%% Confirm your results
% e.g. do ARC proteins consistently fail to co-elute?

[x1,x2,x3] = unique(ds);
tmp = nan(length(x1),size(ravg_corum,2));
for ii = 1:size(tmp,1)
  I = ds==ii;
  tmp(ii,:) = nanmean(ravg_corum(I,:));
end
tmp = ravg_corum;

groups = repmat(1:size(corum,1),size(tmp,1),1);
pp = anovan(tmp(:),groups(:),'display','off');

I = find(sum(~isnan(tmp)) > 0);
[~,I2] = sort(nanmean(tmp(:,I)),'descend');
figure
subplot(2,1,1)
imagesc(tmp(:,I(I2)))
colorbar
title(['p = ' num2str(pp)])
xlabel('Corum complex')
ylabel('Dataset')
subplot(2,1,2)
imagesc(fbad_corum(:,I(I2)))
colorbar
title(['p = ' num2str(pp)])
xlabel('Corum complex')
ylabel('Dataset')

sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/Figure1.eps';
set(gcf,'paperunits','inches','paperposition',[.1 .1 10 5],'units','inches','position',[.1 .1 10 5])
%print(gcf,sf,'-depsc2');



Ndatasets = zeros(size(ravg_corum,2),1);
for ii = 1:length(Ndatasets)
  Ndatasets(ii) = length(unique(ds(~isnan(ravg_corum(:,ii)))));
end
good1 = find(Ndatasets>=1);
tmp = ravg_corum(:,good1);
[~,I] = sort(nanmean(tmp),'ascend');
Ndatasets2 = Ndatasets(good1(I));

figure
subplot(2,1,1)
imagesc(tmp(:,I))
%colorbar
title(['p = ' num2str(pp)])
xlabel('Complex number')
ylabel('Dataset')

subplot(2,1,2),hold on
mm = nanmean(tmp(:,I));
ss = nanstd(tmp(:,I));
for ii = 1:length(I)
  cols = 'b';
  if Ndatasets2(ii)==1
    cols = 'r';
  end
  %plot([ii ii],[prctile(tmp(:,I(ii)),0) prctile(tmp(:,I(ii)),100)],'color',cols,'linestyle',':');
  plot([ii ii],mm(ii)+[-ss(ii) ss(ii)],'color',cols);
end
tmp2 = tmp(:,I);
corumN2 = corumN(I);
scatter(find(Ndatasets2==1),nanmean(tmp2(:,Ndatasets2==1)),corumN2(Ndatasets2==1)'*5,...
  'MarkerFaceColor','r','MarkerEdgeColor','r',...
  'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
scatter(find(Ndatasets2>1),nanmean(tmp2(:,Ndatasets2>1)),corumN2(Ndatasets2>1)'*5,...
  'MarkerFaceColor','b','MarkerEdgeColor','b',...
  'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
%scatter(40,.9,100*5,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
%scatter(40,.75,25*5,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
%scatter(40,.65,2*5,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
%text(60,.9,'100')
%text(55,.75,'25')
%text(53,.65,'2')
plot([700 730],[-.2 -.2],'r')
plot([700 730],[-.1 -.1],'b')
text(735,-.2,'Single dataset')
text(735,-.1,'>1 dataset')
axis([-1 size(tmp,2) -.5 1.01])
ylabel('Average inter-complex pairwise R')
xlabel('Complex number')

sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/Figure2.eps';
set(gcf,'paperunits','inches','paperposition',[.1 .1 12 4],'units','inches','position',[.1 .1 12 4])
%print(gcf,sf,'-depsc2');



%% Now load interactomes from other studies

clear fn
fn{1} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Wan_et_al.tsv';
fn{2} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Havugimana_et_al.csv';
fn{3} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/BioPlex.tsv';
fn{4} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Hein_et_al.tsv';
fn{5} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/HI-II-14.tsv';

interactomes = cell(size(fn));
for ii = 1:length(fn)
  fid = fopen(fn{ii});
  cc = 0;
  interactomes{ii} = cell(100000,2);
  if ~isempty(strfind(fn{ii},'.tsv'))
    separator = '\t';
  elseif ~isempty(strfind(fn{ii},'.csv'))
    separator = ',';
  end
  while ~feof(fid)
    cc = cc+1;
    t1 = strsplit(fgetl(fid),separator);
    interactomes{ii}{cc,1} = t1{1};
    interactomes{ii}{cc,2} = t1{2};
    
  end
  interactomes{ii} = interactomes{ii}(1:cc,:);
end

inter2cor = zeros(8,size(corum,1));
for ii = 1:length(interactomes)
  ii
  for jj = 1:size(interactomes{ii},1)
    I1 = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,1})));
    I2 = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,2})));
    I = intersect(I1,I2);
    
    for kk = 1:length(I)
      inter2cor(ii,I(kk)) = inter2cor(ii,I(kk))+1;
    end
  end
end

% now compare my full datasets
fn1 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D1_Final_Interactions_list_30_precision.csv';
fn2 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D2_Final_Interactions_list_37_precision.csv';
fn3 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D4_Final_Interactions_list_50_precision.csv';
dd{1} = readFinalInteractionList(fn1);
dd{2} = readFinalInteractionList(fn2);
dd{3} = readFinalInteractionList(fn3);
I = dd{1}.data(:,6)>0.5;
dd{1}.data = dd{1}.data(I,:);
dd{1}.text = dd{1}.text(I,:);
I = dd{2}.data(:,6)>0.5;
dd{2}.data = dd{2}.data(I,:);
dd{2}.text = dd{2}.text(I,:);
I = dd{3}.data(:,6)>0.5;
dd{3}.data = dd{3}.data(I,:);
dd{3}.text = dd{3}.text(I,:);

interactomes{6} = dd{1}.text(:,5:6);
interactomes{7} = dd{2}.text(:,5:6);
interactomes{8} = dd{3}.text(:,5:6);

for ii = 6:8
  ii
  for jj = 1:size(dd{ii-5}.text,1)
    protA = dd{ii-5}.text{jj,2};
    protB = dd{ii-5}.text{jj,3};
    I1 = find(~cellfun('isempty',strfind(corum(:,5),protA)));
    I2 = find(~cellfun('isempty',strfind(corum(:,5),protB)));
    I = intersect(I1,I2);
    
    for kk = 1:length(I)
      inter2cor(ii,I(kk)) = inter2cor(ii,I(kk))+1;
    end
  end
end



%% What CORUM complexes do other studies find? (nice figure)

fns = {'Wan (CE-MS)', 'Havug. (CE-MS)', 'Bioplex (AP-MS)','Hein (AP-MS)','Rolland (Y2H)', ...
  'Nick, apoptosis' 'Nick, HeLa' 'Anders'};
fns_reduced = {'CEMS' 'CEMS' 'APMS' 'APMS' 'Y2H' 'N.A.' 'N.H.' 'And.'};

for mm = 1:2
  
  figure
  
  if mm==1
    [~,II] = sort(nansum(inter2cor(6:8,:)>0),'descend');
  elseif mm==2
    [~,II] = sort(nansum(inter2cor(6:8,:)>0),'descend');
    good1 = find(Ndatasets>2);
    II = II(ismember(II,good1));
  end
  
  subplot(5,2,1:4)
  %imagesc(double(inter2cor(1:8,II)>0))
  imagesc(log10(inter2cor(1:8,II)))
  colorbar
  %y=get(colorbar,'YTick');
  y = [0 1 2 3];
  colorbar('ytick',y,'YTickLabel',10.^y);
  set(gca,'ytick',1:8,'yticklabel',fns_reduced)
  xlabel('CORUM complex number')
  if mm==1
    title('Number of interactions predicted, all CORUM complexes')
  elseif mm==2
    title('Number of interactions predicted, detected CORUM complexes')
  end
  
  % add random rows to inter2cor
  for ii = 9:500
    inter2cor(ii,II) = inter2cor(mod(ii,3)+6,randsample(II,length(II)));
    inter2cor(ii,II) = inter2cor(mod(ii+1,3)+6,randsample(II,length(II)));
    inter2cor(ii,II) = inter2cor(mod(ii+2,3)+6,randsample(II,length(II)));
  end
  
  [R,Rp] = corr(inter2cor(:,II)','type','pearson');
  for ii = 1:size(R,1); R(ii,ii) = nan;end
  subplot(5,2,[5 7])
  imagesc(R(1:11,1:11))
  caxis([nanmin(R(:)) nanmax(R(:))])
  %set(gca,'xtick',1:8,'xticklabel',fns_reduced)
  colorbar
  axis square
  
  subplot(5,2,9),hold on
  xx = [.95 1.05 1.95 2.05 3];
  yR = nanmean(R(1:5,6:8),2);
  yRchance = nanmean(R(1:5,9:end),2);
  yRchance_sd = nanstd(R(1:5,9:end),[],2);
  scatter(xx,yR,20,'b','filled')
  scatter(xx,yRchance,20,'r','filled')
  for ii = 1:length(xx)
    plot(xx([ii ii]),yRchance(ii) + yRchance_sd*[-1 1],'-r')
  end
  ylabel('Avg. R')
  set(gca,'xtick',1:3,'xticklabel',{'Co-Elution' 'AP-MS' 'Y2H'})
  xlim([.5 3.5])
  %title('Correlation')
  
  E = nan(11,11);
  Ep = nan(11,11);
  Ntotal = length(II);
  for ii = 1:8+3
    if ii>11 % random chance
      com1 = find(inter2cor(ii-3,randsample(II,length(II)))>0);
    else
      com1 = find(inter2cor(ii,II)>0);
    end
    frac1 = [length(com1) Ntotal];
    for jj = 1:8+3
      if ii==jj; continue; end
      if jj>11
        com2 = find(inter2cor(jj-3,randsample(II,length(II)))>0);
      else
        com2 = find(inter2cor(jj,II)>0);
      end
      
      I = intersect(com1,com2);
      
      frac2 = [length(I) length(com2)];
      Ep(ii,jj) = mychiproptest(frac1(1),frac1(2),frac2(1),frac2(2));
      E(ii,jj) = (frac2(1) / frac2(2)) / (frac1(1) / frac1(2));
    end
  end
  subplot(5,2,[6 8])
  imagesc(E)
  caxis([0 nanmax(E(:))])
  set(gca,'xtick',1:8,'ytick',1:8);%'xticklabel',fns_reduced,'yticklabel',fns_reduced)
  colorbar
  axis square
  %title('Enrichment')
  
  subplot(5,2,10),hold on
  yE = nanmean(E(6:8,1:5),1);
  yEchance = nanmean(E(1:5,9:end),2);
  yEchance_sd = nanstd(E(1:5,9:end),[],2);
  scatter(xx,yE,20,'b','filled')
  scatter(xx,yEchance,20,'r','filled')
  for ii = 1:length(xx)
    plot(xx([ii ii]),yEchance(ii) + yEchance_sd*[-1 1],'-r')
  end
  ylabel('Avg. Enrichment')
  set(gca,'xtick',1:3,'xticklabel',{'Co-Elution' 'AP-MS' 'Y2H'})
  xlim([.5 3.5])
  
  set(gcf,'units','normalized','position',[.3 .1 .5 .75])
  
end



%% What CORUM complexes do other studies find? (try presenting it another way)

fns = {'Wan (CE-MS)', 'Havug. (CE-MS)', 'Bioplex (AP-MS)','Hein (AP-MS)','Rolland (Y2H)', ...
  'Nick, apoptosis' 'Nick, HeLa' 'Anders'};
fns_reduced = {'CEMS' 'CEMS' 'APMS' 'APMS' 'Y2H' 'Nick1' 'Nick2' 'Anders'};

I = [6 7 8 1 2 3 4 5];
inter2cor2 = inter2cor(I,:);
fns = fns(I);
fns_reduced = fns_reduced(I);

[~,II] = sort(nansum(inter2cor2(1:3,:)),'descend');
good1 = find(Ndatasets>2);
II = II(ismember(II,good1));

% enrichment analysis
% N = all corum
% A = corum complexes that we found
% B = corum complexes found by a dataset
% AB = complexes found by us and the dataset
% compare A/N to AB/B
iterMax = 100;
E = nan(1,size(inter2cor2,1));
Ep = nan(1,size(inter2cor2,1));
E_chance = nan(iterMax,size(inter2cor2,1));
RR = nan(size(E));
A = find(nansum(inter2cor2(1:3,II))>0); % our CORUM complexes
for ii = 1:size(inter2cor2,1)
  B = find(inter2cor2(ii,II)>0); % CORUM complexes in this dataset
  AB = intersect(A,B);
  E(ii) = (length(AB)/length(B)) / (length(A)/length(II));
  Ep(ii) = mychiproptest(length(AB),length(B),length(A),length(II));
  
  % measure chance enrichment, too
  for iter = 1:iterMax
    tmp = randsample(length(II),length(B));
    AB = intersect(A,tmp);
    E_chance(iter,ii) = (length(AB)/length(tmp)) / (length(A)/length(II));
  end
  
  % correlation
  RR(ii) = corr(nansum(inter2cor2(1:3,:))',inter2cor2(ii,:)','type','Pearson');
end

figure
subplot(2,1,1)
imagesc(log10(inter2cor2(1:8,II)))
colorbar
%y=get(colorbar,'YTick');
y = [0 1 2 3];
colorbar('ytick',y,'YTickLabel',10.^y);
set(gca,'ytick',1:8,'yticklabel',fns,'fontsize',8)
xlabel('CORUM complex number')
ylabel('Dataset')
title('Number of predicted interactions')

subplot(2,2,4),hold on
scatter(1:length(E),E,50,'k','filled')
yy = nanmean(E_chance);
ys = nanstd(E_chance);
myplotpatch(0:length(E)+1,yy([1 1:end end]),ys([1 1:end end]),[1 .9 .9])
xlim([.5 length(E)+.5])
ylabel('Enrichment value')
set(gca,'xtick',1:8,'xticklabel',fns,'XTickLabelRotation',25,'fontsize',8)
title('Enrichment')

subplot(2,2,3),hold on
scatter(1:length(RR),RR,50,'k','filled')
xlim([.5 length(E)+.5])
set(gca,'xtick',1:8,'xticklabel',fns,'XTickLabelRotation',25,'fontsize',8)
ylabel('Pearson R')
title('Correlation')

set(gcf,'units','normalized','position',[.1 .1 .8 .4],...
  'paperunits','normalized','paperposition',[.1 .1 .8 .4])
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/compareDatasets';
%print(gcf,sf,'-depsc2');
print(gcf,sf,'-dpng');


%% Focus on a few complexes

% show the connection matrices for ribosome from the different studies...
% We, and other CE-MS, finda  lot of ribosome interactions.
% The matrix should be full for us, and similarly for other CE-MS.
% Hopefully similar patterns too.
% It shoud be emptier for AP-MS and Y2H.


[Ispliceosome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'Spliceosome')));
[Iribosome,x1] = find(~cellfun('isempty',strfind(corum,'55S ribosome')));
[Iproteasome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'20S proteasome')));
[Isignalosome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'COP9 signalosome complex (G')));

%I2 = {Ispliceosome(1) Iribosome(1) Iproteasome(1) Isignalosome(1)};

I = [6 7 8 1 2 3 4 5];
interactomes2 = interactomes(I);

for kk = 1
  I2 = II(kk);
  Prots_ribosome = strsplit(corum{I2,5},',');
  C_ribosome = cell(1,8);
  for ii = 1:length(interactomes2)
    C_ribosome{ii} = nan(length(Prots_ribosome),length(Prots_ribosome));
    
    ia = cell(length(Prots_ribosome),1);
    for jj = 1:length(Prots_ribosome)
      [ia{jj},ib] = find(ismember(interactomes2{ii},Prots_ribosome{jj}));
    end
    for jj = 1:length(Prots_ribosome)
      for kk = 1:length(Prots_ribosome)
        if jj==kk; continue; end
        C_ribosome{ii}(jj,kk) = ~isempty(intersect(ia{jj},ia{kk}));
      end
    end
  end
  
  complexName = corum{I2,2};
  
  figure
  for ii = 1:8
    subplot(2,4,ii)
    imagesc(C_ribosome{ii})
    axis square
    caxis([0 1])
    title(fns{ii},'fontsize',6)
    set(gca,'fontsize',5)
    ylabel('Protein number','fontsize',5)
    xlabel('Protein number','fontsize',5)
  end
  set(gcf,'units','normalized','position',[.1 .1 .8 .4],...
    'paperunits','normalized','paperposition',[.1 .1 .8 .4])
  sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/compareDatasets2';
  %print(gcf,sf,'-depsc2');
  print(gcf,sf,'-dpng');
  
  disp(corum{I2,2})
end