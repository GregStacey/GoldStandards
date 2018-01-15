
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


%% Load data

% 1. Individual replicates

%datasets = {'Anders' 'Nick_HeLa' 'Nick_Ap' 'Nick_tiss' 'Jenny' 'Craig' 'Mike' 'Mike_{10}'};
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
clear X Dist TP_Matrix binary_interaction_list A inverse_self possibleInts

Ndatasets = zeros(size(ravg_corum,2),1);
for ii = 1:length(Ndatasets)
  Ndatasets(ii) = length(unique(ds(~isnan(ravg_corum(:,ii)))));
end
good1 = find(Ndatasets>=1);
tmp = ravg_corum(:,good1);
[~,I] = sort(nanmean(tmp),'ascend');
Ndatasets2 = Ndatasets(good1(I));


% 2. Full interaction lists from our lab
interactome_description = {'Us' 'Us' 'Us' ...
  'ce' 'ce' 'ce' ...
  'ap-ms' 'ap-ms' 'ap-ms' ...
  'y2h' 'y2h' 'y2h'};
Ilargescale = [1 1 1 1 1 1 1 1 0 1 0 1 0 1]==1;

fn1 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D1_Final_Interactions_list_30_precision.csv';
fn2 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D2_Final_Interactions_list_37_precision.csv';
fn3 = '/Users/Mercy/Academics/Foster/Manuscripts/PCPSILAC_Methods/Data/D4_Final_Interactions_list_50_precision.csv';
clear dd
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

interactomes = cell(12,1);
interactomes{1} = dd{1}.text(:,5:6);
interactomes{2} = dd{2}.text(:,5:6);
interactomes{3} = dd{3}.text(:,5:6);


% 3. Published interactions lists from other labs

clear fn
fn{1} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Wan_et_al.tsv';            % co-elution
fn{2} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Havugimana_et_al.csv';     % co-elution
fn{3} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Kirkwood_et_al.csv';     % co-elution
fn{4} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/BioPlex.tsv';              % ap-ms
fn{5} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Hein_et_al.tsv';           % ap-ms
fn{6} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Ewing_et_al.csv';          % ap-ms
fn{7} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/HI-II-14.tsv';             % y2h
fn{8} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Wang_et_al.csv';           % y2h
fn{9} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Rual_et_al.csv';          % y2h
%fn{6} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Varjosalo_et_al.csv';      % ap-ms
%fn{8} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Stelzl_et_al.csv';         % y2h
%fn{8} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Bandyopadhyay_et_al.csv'; % y2h

for ii = 1:length(fn)
  ii+3
  fid = fopen(fn{ii});
  cc = 0;
  interactomes{ii+3} = cell(100000,2);
  if ~isempty(strfind(fn{ii},'.tsv'))
    separator = '\t';
  elseif ~isempty(strfind(fn{ii},'.csv'))
    separator = ',';
  end
  while ~feof(fid)
    cc = cc+1;
    t1 = strsplit(fgetl(fid),separator);
    interactomes{ii+3}{cc,1} = t1{1};
    interactomes{ii+3}{cc,2} = t1{2};
    
  end
  interactomes{ii+3} = interactomes{ii+3}(1:cc,:);
end


% 4. Match interactomes to corum

% make empty connection matrices
connMat = cell(length(interactomes),size(corumN,1));
for ii = 1:length(interactomes)
  for jj = 1:size(corumN,1)
    connMat{ii,jj} = zeros(corumN(jj));
  end
end

for ii = 1:length(interactomes)
  ii
  for jj = 1:size(interactomes{ii},1)
    % no self-interactions
    if strcmp(interactomes{ii}{jj,1},interactomes{ii}{jj,2});
      continue;
    end
    I1 = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,1})));
    I2 = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,2})));
    I = intersect(I1,I2);
    
    for kk = 1:length(I)
      protNames = strsplit(corum{I(kk),5},',');
      I3 = find(ismember(protNames,interactomes{ii}{jj,1}));
      I4 = find(ismember(protNames,interactomes{ii}{jj,2}));
      connMat{ii,I(kk)}(I3,I4) = 1;
      connMat{ii,I(kk)}(I4,I3) = 1;
    end
  end
end

inter2cor = zeros(size(interactomes,1),size(corum,1));
for ii = 1:length(interactomes)
  for jj = 1:size(corum,1)
    inter2cor(ii,jj) = nansum(connMat{ii,jj}(:)) / 2;
  end
end




%% What CORUM complexes do other studies find?

[~,II] = sort(nansum(inter2cor(1:3,:)),'descend');
good1 = find(Ndatasets>2);
II = II(ismember(II,good1));
inter2cor2 = inter2cor(:,II);

corCov = nan(size(inter2cor2));
for ii = 1:length(interactomes)
  corCov(ii,:) = inter2cor2(ii,:) ./ corumN(II)' ./ (corumN(II)-1)' * 2;
end

% enrichment analysis
% N = all corum
% A = corum complexes that we found
% B = corum complexes found by someone else
% AB = complexes found by us and them
% compare A/N to AB/B
E = nan(length(interactomes),length(interactomes));
Ep = nan(length(interactomes),length(interactomes));
NN = nan(length(interactomes),length(interactomes));
Ntotal = length(II);
for ii = 1:length(interactomes)
  com1 = find(inter2cor(ii,II)>0);
  frac1 = [length(com1) Ntotal];
  for jj = 1:length(interactomes)
    if ii==jj; continue; end
    com2 = find(inter2cor(jj,II)>0);
    I = intersect(com1,com2);
    
    if length(I)>1
      frac2 = [length(I) length(com2)];
      Ep(ii,jj) = mychiproptest(frac1(1),frac1(2),frac2(1),frac2(2));
      E(ii,jj) = (frac2(1) / frac2(2)) / (frac1(1) / frac1(2));
      NN(ii,jj) = length(I);
    end
  end
end

RR2 = nan(length(interactomes));
RRp = nan(length(interactomes));
RR3 = nan(length(interactomes));
for ii = 1:length(interactomes)
  for jj = 1:length(interactomes)
    if ii==jj; continue; end
    I = inter2cor2(ii,:)>0 & inter2cor2(jj,:)>0;
    if sum(I)>1
      [RR2(ii,jj),RRp(ii,jj)] = corr(inter2cor2(ii,I)',inter2cor2(jj,I)','type','spearman');
      RR3(ii,jj) = corr(corCov(ii,I)',corCov(jj,I)','type','pearson');
    end
  end
end



figure,colormap bone
subplot(3,3,1:3)
imagesc(log10(inter2cor2))
colorbar
y = [0 1 2 3];
colorbar('ytick',y,'YTickLabel',10.^y);
set(gca,'ytick',1:length(interactomes),'yticklabel',interactome_description,'fontsize',11)
xlabel('CORUM complex number')
ylabel('Dataset')
title('Number of interactions per complex')

subplot(3,3,4:6)
imagesc(corCov)
colorbar
y = [0 1 2 3];
colorbar
%colorbar('ytick',y,'YTickLabel',10.^y);
set(gca,'ytick',1:length(interactomes),'yticklabel',interactome_description,'fontsize',11)
xlabel('CORUM complex number')
ylabel('Dataset')
title('Coverage per complex')

subplot(3,3,7)
imagesc(RR3)
colorbar
set(gca,'xtick',1:length(interactomes),'ytick',1:length(interactomes),...
  'xticklabel',interactome_description,'yticklabel',interactome_description)
axis square
title('Correlation')
caxis([0 1])

subplot(3,3,8)
tmpE = (Ep);
tmpE(tmpE==0) = nan;
pp = sort(tmpE(:));
pp(isnan(pp)) = [];
tmp = (1:length(pp)) / length(pp) *.05;
J = find(pp'-tmp<0,1,'last');
xcut = 0;
if ~isempty(J)
  xcut = pp(J);
end
imagesc(Ep<=xcut)
colorbar
set(gca,'xtick',1:length(interactomes),'ytick',1:length(interactomes),...
  'xticklabel',interactome_description,'yticklabel',interactome_description)
axis square
title('Significant enrichment')

subplot(3,3,9)
loglog(1,1,'k'),hold on
loglog(1,1,'b')
loglog(1,1,'r')
loglog(1,1,'m')
cols = {'k' 'b' 'r' 'm'};
ss = {'Us' 'ce' 'ap-ms' 'y2h'};
for ii = 1:length(interactomes)
  yy = sort(inter2cor(ii,:),'descend');
  %yy2 = sort(inter2cor(ii,:),'descend') / length(interactomes{ii});
  I = find(ismember(ss,interactome_description{ii}));
  loglog(yy,cols{I})
end
grid on
axis([0 1000 0 3000])
title('Size of complexes')
xlabel('Complex number')
ylabel('Predicted complex size')
legend('us','ce','ap-ms','y2h','location','northeast')

set(gcf,'units','normalized','position',[0 0 .7 1])
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/Figure3';
%print(gcf,sf,'-depsc2');



%% relegated figures

% figure
% ss = {'Us' 'ce' 'ap-ms' 'y2h'};
% for ii = 1:4
%   for jj = 1:4
%     I1 = find(ismember(interactome_description,ss{ii}));
%     I2 = find(ismember(interactome_description,ss{jj}));
%
%
%     subplot(4,4, (ii-1)*4 + jj)
%     scatter(nanmean(corCov(I1,II)),nanmean(corCov(I2,II)),corumN(II)'*5,...
%       'MarkerFaceColor','b','MarkerEdgeColor','b',...
%       'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2)
%     axis([0 1 0 1])
%   end
% end

% RR3 = nan(4,4);
% ss = {'Us' 'ce' 'ap-ms' 'y2h'};
% for ii = 1:4
%   for jj = 1:4
%     %if ii ==jj; continue; end
%     I1 = ismember(interactome_description,ss{ii});
%     I2 = ismember(interactome_description,ss{jj});
%
%     A = nanmean(inter2cor2(I1,:));
%     B = nanmean(inter2cor2(I2,:));
%     RR3(ii,jj) = corr(A',B','type','pearson');
%   end
% end



%% Focus on a individual connection matricies

% show the connection matrices for ribosome from the different studies...
% We, and other CE-MS, finda  lot of ribosome interactions.
% The matrix should be full for us, and similarly for other CE-MS.
% Hopefully similar patterns too.
% It shoud be emptier for AP-MS and Y2H.


%[Ispliceosome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'Spliceosome')));
%[Iribosome,x1] = find(~cellfun('isempty',strfind(corum,'55S ribosome')));
%[Iproteasome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'20S proteasome')));
%[Isignalosome,x1] = find(~cellfun('isempty',strfind(corum(:,2),'COP9 signalosome complex (G')));

%I2 = {Ispliceosome(1) Iribosome(1) Iproteasome(1) Isignalosome(1)};

%I = [6 7 8 1 2 3 4 5];
%interactomes2 = interactomes(I);

[~,II] = sort(nansum(inter2cor(1:3,:)),'descend');

subplots = [1 2 3 4 5 6 7 8 9 10 11 12];

for kk = 1%:100
  
  figure
  for ii = 1:length(interactomes)
    cc = subplots(ii);
    subplot(4,3,cc)
    imagesc(connMat{ii,II(kk)})
    axis square
    caxis([0 1])
    title(interactome_description{ii},'fontsize',6)
    set(gca,'fontsize',8)
  end
  set(gcf,'units','inches','position',[0 0 12 12],...
    'paperunits','inches','paperposition',[0 0 12 12])
  sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/compareDatasets2';
  %print(gcf,sf,'-depsc2');
  
  disp(corum{II(kk),2})
  %pause
  %close all
end

% do connection matrices look similar?
[~,II] = sort(nansum(inter2cor(1:3,:)),'descend');
good1 = find(Ndatasets>1);
II = II(ismember(II,good1));
RR4 = zeros(length(interactomes));
for ii = 1:length(interactomes)
  for jj = 1:length(interactomes)
    if ii==jj; continue; end
    rr = zeros(length(II),1);
    for kk = 1:length(II)
      A = find(connMat{ii,II(kk)}(:)==1);
      B = find(connMat{jj,II(kk)}(:)==1);
      nn = length(intersect(A(:),B(:)));
      rr(kk) = nn / (length(A) + length(B) + nn);
    end
    
    RR4(ii,jj) = nanmean(rr);
  end
end

II = II(1:50);
RR5 = zeros(length(interactomes));
for ii = 1:length(interactomes)
  for jj = 1:length(interactomes)
    if ii==jj; continue; end
    rr = zeros(length(II),1);
    for kk = 1:length(II)
      A = find(connMat{ii,II(kk)}(:)==1);
      B = find(connMat{jj,II(kk)}(:)==1);
      nn = length(intersect(A(:),B(:)));
      rr(kk) = nn / (length(A) + length(B) + nn);
    end
    
    RR5(ii,jj) = nanmean(rr);
  end
end

figure
subplot(1,2,1)
imagesc(RR4),colorbar
caxis([0 .25])
axis square
set(gca,'xtick',1:length(interactomes),'ytick',1:length(interactomes),...
  'xticklabel',interactome_description,'yticklabel',interactome_description)
title('550 detectable complexes')
subplot(1,2,2)
imagesc(RR5),colorbar
caxis([0 .25])
axis square
set(gca,'xtick',1:length(interactomes),'ytick',1:length(interactomes),...
  'xticklabel',interactome_description,'yticklabel',interactome_description)
title('50 best complexes')

set(gcf,'units','inches','position',[0 0 12 5],...
  'paperunits','inches','paperposition',[0 0 12 5])
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/compareDatasets3';
%print(gcf,sf,'-depsc2');



%% What's special about the 50?

% get annotation for human proteins
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/goa_human.gaf';
A = goannotread(sf,'Fields',{'DB_Object_ID','GOid'});

% all GO terms
GO = geneont('live',true);

% collate all corum proteins
[~,II] = sort(nansum(inter2cor(1:3,:)),'descend');
good1 = find(Ndatasets>2);
II = II(ismember(II,good1));
II = II(1:50);
proteins = cell(10^5,1);
the50 = nan(10^5,2); % 1059-by-2, indices of the50 proteins in 'proteins'
kk1 = 0;
kk2 = 0;
for ii = 1:size(corum,1)
  tmp = strsplit(corum{ii,5},',');
  I = kk1+1 : kk1+length(tmp);
  proteins(I) = tmp;
  kk1 = kk1+length(tmp);
  
  if ismember(ii,II)
    I2 = kk2+1 : kk2+length(tmp);
    the50(I2,1) = I;
    the50(I2,2) = ii;
    kk2 = kk2+length(tmp);
  end
end
proteins = proteins(1:kk1);
the50 = the50(1:kk2,:);

% for efficiency, make a map of Protein_ID --> GOid
protID = cell(length(A),1);
goID = cell(length(A),1);
for ii = 1:length(A)
  protID{ii} = A(ii).DB_Object_ID;
  goID{ii} = A(ii).GOid;
end

% reduce Annotations and Proteins to just intersecting members
[~,I2a,I2p] = intersect(protID,proteins);
A2 = A(I2a);
proteins2 = proteins(I2p);

id2goMap = containers.Map;
for ii = 1:numel(A2)
  key = A2(ii).DB_Object_ID;
  if isKey(id2goMap,key)
    id2goMap(key) = [id2goMap(key) A2(ii).GOid];
  else
    id2goMap(key) = A2(ii).GOid;
  end
end


% count the number of annotations for
%   - every protein
%   - proteins in the 50
mm = GO.Terms(end).id;           % gets the last term id
proteinscount = zeros(mm,1);     % a vector of GO term counts for all corum proteins.
the50count = zeros(mm,1);        % a vector of GO term counts for the 50.
the50count_ind = zeros(mm,50);
for ii = 1:length(proteins2)
  ii
  if isKey(id2goMap,proteins2{ii})
    goid = getrelatives(GO,id2goMap(proteins2{ii}));
    proteinscount(goid) = proteinscount(goid) + 1;
    
    I = the50(:,1)==I2p(ii);
    if sum(I)>0
      the50count(goid) = the50count(goid) + 1;
      
      corumI = II == the50(I,2);
      the50count_ind(goid,corumI) = the50count_ind(goid,corumI) + 1;
    end
  end
end

% terms enriched in The 50
I = sum(the50count_ind,2)>0;

pp = zeros(51,sum(I));
pp(1,:) = hygepdf(the50count(I),max(proteinscount(I)),...
  max(the50count(I)),proteinscount(I));
for ii = 1:50
  % terms enriched in any single complex of The 50
  pp(ii+1,:) = hygepdf(the50count_ind(I,ii),max(proteinscount(I)),max(the50count_ind(I,ii)),proteinscount(I));
end

% B-H correction (only on pvalues != 1)
Ncomp = sum(pp(:)~=1);
[~,I2] = sort(pp(:));
bhcomp = (1:Ncomp) / Ncomp * .05;
I3 = find(bhcomp'-pp(I2(1:Ncomp))>0,1,'last');
[ia,ib] = ind2sub(size(pp),I2(1:I3));


% Find annotation terms that are sig. for The 50, but not individual complexes

sigTerms = pp(1,:) < pp(I2(I3));
nn = sum(pp(2:end,:) < pp(I2(I3)));

tmp = find(sigTerms & nn<2);
I = find(sum(the50count_ind,2)>0);

disp('INTERESTING GO TERMS')
for ii = 1:length(tmp)-5
  disp(GO.Terms(I(tmp(ii))).name)
end


%%

[~,II] = sort(nansum(inter2cor(1:3,:)),'descend');
good1 = find(Ndatasets>2);
II = II(ismember(II,good1));
II = II(1:50);


% Purification method?
S1 = corum(:,7);

% collate "purification IDs"
purID = zeros(size(S1));
purDesc = cell(size(S1));
ID50 = zeros(size(S1));
kk = 0;
for ii = 1:size(purID)
  s1 = strsplit(S1{ii},'|');
  for jj = 1:length(s1)
    kk = kk+1;
    I = find(s1{jj} == ':');
    purID(kk) = str2double(s1{jj}(I+1:I+4));
    
    purDesc{kk} = s1{jj}(I+7 : end);
    
    ID50(kk) = 0;
    if ismember(ii,II)
      ID50(kk) = 1;
    end
  end
%   I = find(S1{ii} == ':');
%   for jj = 1:length(I)
%     kk = kk+1;
%     purID(kk) = str2double(S1{ii}(I+1:I+4));
%     
%     ID50(kk) = 0;
%     if ismember(ii,II)
%       ID50(kk) = 1;
%     end
%   end
end

[unqID,I2] = unique(purID);
unqDesc = purDesc(I2);
pp = zeros(size(unqID));
EE = zeros(size(pp));
for ii = 1:length(unqID)
  frac1 = [sum(purID(ID50==1)==unqID(ii)) length(purID(ID50==1))];
  frac2 = [sum(purID==unqID(ii)) length(purID)];
  EE(ii) = (frac1(1)/frac1(2)) / (frac2(1)/frac2(2));
  pp(ii) = mychiproptest(frac1(1),frac1(2),frac2(1),frac2(2));
end

bhcomp = (1:length(pp)) / length(pp) * .01;

[~,I] = sort(pp);
nsig = find(bhcomp'-sort(pp)>0,1,'last');

ss = sprintf('\nThe 50 are enriched for ... (purification method)');
disp(ss)
for ii = 1:nsig
  ss1 = [unqDesc{I(ii)} ',\t enrichment = ' num2str(EE(I(ii))) ' (p = ' num2str(pp(I(ii))) ')'];
  ss2 = sprintf(ss1);
  disp(ss2)
end




% species
[unqSpecies,Is,Is2] = unique(corum(:,4));
pp2 = zeros(size(Is));
EE2 = zeros(size(Is));
for ii = 1:length(unqSpecies)
  frac1 = [sum(Is2(II)==ii) length(II)];
  frac2 = [sum(Is2==ii) length(Is2)];
  EE2(ii) = (frac1(1)/frac1(2)) / (frac2(1)/frac2(2));
  pp2(ii) = mychiproptest(frac1(1),frac1(2),frac2(1),frac2(2));
end

I = find(pp2<.05/length(pp2));
ss = sprintf('\nThe 50 are enriched for ... (species)');
disp(ss)
for ii = 1:length(I)
  ss1 = [unqSpecies{I(ii)} ',\t enrichment = ' num2str(EE2(I(ii))) ' (p = ' num2str(pp2(I(ii))) ')'];
  ss2 = sprintf(ss1);
  disp(ss2)
end
