
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


%% Read in and translate Affinity Benchmark

% write PDB ids to translate to uniprot
% fn = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/data/affinity_benchmark1.tsv';
% fid = fopen(fn);
% fgetl(fid);
% head = strsplit(fgetl(fid),'\t','collapsedelimiters',0);
% fnout = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/data/AB1_pdb.csv';
% fidout = fopen(fnout,'w');
% while ~feof(fid)
%   tt = fgetl(fid);
%   tt = strrep(tt,'"','');
%   t1 = strsplit(tt,'\t','collapsedelimiters',0);
%
%   protA = strsplit(t1{3},'_');
%   protA = protA{1};
%   protB = strsplit(t1{5},'_');
%   protB = protB{1};
%
%   fprintf(fidout,'%s,%s\n',protA,protB);
% end

% read Affinity Benchmark
fn = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/data/affinity_benchmark1.tsv';
fid = fopen(fn);
fgetl(fid);
head = strsplit(fgetl(fid),'\t','collapsedelimiters',0);
AB.data = nan(10^4,3);
AB.text = cell(10^4,5);
cc = 0;
while ~feof(fid)
  tt = fgetl(fid);
  tt = strrep(tt,'"','');
  t1 = strsplit(tt,'\t','collapsedelimiters',0);
  
  protA = strsplit(t1{3},'_');
  protA = protA{1};
  protB = strsplit(t1{5},'_');
  protB = protB{1};
  
  cc = cc+1;
  AB.data(cc,1) = str2double(t1{8});  % Kd
  AB.data(cc,2) = str2double(t1{9}); % dG
  AB.data(cc,3) = str2double(t1{10}); % I-RMSD
  AB.data(cc,4) = str2double(t1{11}); % DASA
  AB.text{cc,3} = protA;
  AB.text{cc,4} = protB;
  AB.text{cc,5} = t1{2};
end
AB.data = AB.data(1:cc,:);
AB.text = AB.text(1:cc,:);

% map from PDB to uniprot
fnf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/data/AB1_uniprot.tsv';
fid = fopen(fnf);
head = strsplit(fgetl(fid),'\t');
cc = 0;
while ~feof(fid)
  t2 = strsplit(fgetl(fid),'\t');
  poss_pdb = strsplit(t2{8},';');
  uniprot = t2{1};
  
  Iab = find(ismember(AB.text(:,3),poss_pdb));
  for jj = 1:length(Iab)
    AB.text{Iab(jj),1} = uniprot;
  end
  
  Iab = find(ismember(AB.text(:,4),poss_pdb));
  for jj = 1:length(Iab)
    AB.text{Iab(jj),1} = uniprot;
  end
end



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
fn{7} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/HI-II-14.tsv';             % y2h, Rolland, ..., Vidal 2014
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

for ii = 1:length(interactomes)
  for jj = 1:size(interactomes{ii},1)
    protA = interactomes{ii}{jj,1};
    protB = interactomes{ii}{jj,2};
    tmp = sort({protA protB});
    interactomes{ii}{jj,3} = [tmp{1} '_' tmp{2}];
  end
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


% 5. Find the "enrichment value" for each CORUM complex
EE_complex = nan(length(interactomes),size(corum,1));
EEp_complex = nan(length(interactomes),size(corum,1));
for ii = 1:length(interactomes)
  tmp = interactomes{ii}(:,1:2);
  unqProts = unique(tmp);
  Npotential = length(unqProts) * (length(unqProts)-1) / 2;
  frac_background = [length(interactomes{ii})  Npotential];
  for jj = 1:size(corum,1)
    if mod(jj,250)==1
      disp(['Interactome ' num2str(ii) ', Complex ' num2str(jj) ' of ' num2str(size(corum,1))])
    end
    protNames = strsplit(corum{jj,5},',');
    Nprots_in_interactome = length(intersect(protNames,unqProts));
    if Nprots_in_interactome<2; continue; end
    Nints_complex = sum(connMat{ii,jj}(:)/2);
    Npotential_complex = Nprots_in_interactome * (Nprots_in_interactome-1) / 2;
    frac_complex = [Nints_complex Npotential_complex];
    EEp_complex(ii,jj) = mychiproptest2(frac_complex(1),frac_complex(2),frac_background(1),frac_background(2));
    EE_complex(ii,jj) = (frac_complex(1) / frac_complex(2)) / (frac_background(1) / frac_background(2));
  end
end
EEp_complex(EEp_complex==0) = 1.4822e-323;


% 6. (just for curiosity) What proportion of each interactome is in corum?
incorumN = nan(size(interactomes));
for ii = 1:length(interactomes)
  ii
  incorum = zeros(length(interactomes{ii}),1);
  for jj = 1:length(interactomes{ii})
    ia = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,1})));
    ib = find(~cellfun('isempty',strfind(corum(:,5),interactomes{ii}{jj,2})));
    incorum(jj) = length(intersect(ia,ib));
  end
  incorumN(ii) = sum(incorum) / length(incorum);
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




%% Using complex-enrichment, what complexes are CONSISTENTLY OVER-represented?

disp(' ')

% use subset of corum that's actually seen by some interactomes
%Igood = find(sum(inter2cor>0) > 0);
Igood = find(sum(inter2cor>0) > 0 & corumN'>4);

% find p-value threshold via BH
pvals = EEp_complex(:,Igood);
pvals = pvals(:);
pvals(isnan(pvals)) = [];
pvals = sort(pvals,'ascend');
comp = (1:length(pvals)) / length(pvals) * .01;
ia = find(pvals-comp'>0,1,'first');
pthresh = pvals(ia);

% find consistency
%types = {'Us' 'ce' 'ap-ms' 'y2h'};
%Iconsistent = cell(size(types));
%for ii = 1:length(types)
%  I = ismember(interactome_description,types{ii});
%  Iconsistent{ii} = find(sum(EEp_complex(I,:) < pthresh) >= 2);
%end
Iconsistent{1} = find(sum(EEp_complex(1:3,:) < pthresh) == 3);
Iconsistent{2} = find(sum(EEp_complex(4:5,:) < pthresh) == 2);
Iconsistent{3} = find(sum(EEp_complex(7:9,:) < pthresh) == 3);
Iconsistent{4} = find(sum(EEp_complex(10:12,:) < pthresh) == 3);

% Found in every dataset?
I = mintersect(Iconsistent{1},Iconsistent{2},Iconsistent{3},Iconsistent{4});
nn = 0;
for ii = 1:length(I)
  nc = corumN(I(ii));
  nn = nn+nc*(nc-1)/2;
end
disp([num2str(length(I)) ' complexes found in every dataset (' num2str(nn) ' pairwise interactions):'])
ee = nanmean(EE_complex(:,I));
[~,Isort] = sort(ee,'descend');
I = I(Isort);
ee = round(ee(Isort));
for ii = 1:length(I)
  A = zeros(size(connMat{1,I(ii)}));
  nn = size(A,1);
  for jj = 1:length(interactomes)
    A = A+connMat{jj,I(ii)};
  end
  avgDens = nansum(A(:))/2 / (nn * (nn-1) /2) / length(interactomes);
  complexName = corum{I(ii),2};
  complexName = complexName(1:min([15 length(complexName)]));
  ss = sprintf('%s, %dx enriched, %1.2f avg density,  %d members,  %s',complexName,ee(ii),avgDens,corumN(I(ii)),corum{I(ii),4});
  disp(ss)
end
disp(' ')

% In every co-elution?
%I = intersect(Iconsistent{1},Iconsistent{2});
%I = find(sum(EEp_complex(1:6,:)<pthresh)>=5);
I = find(sum(EEp_complex(1:5,:)<pthresh)>=5 & corumN'>0);
nn = 0;
for ii = 1:length(I)
  nc = corumN(I(ii));
  nn = nn+nc*(nc-1)/2;
end
disp([num2str(length(I)) ' complexes found in every good co-elution dataset (' num2str(nn) ' pairwise interactions):'])
ee = nanmean(EE_complex(1:6,I));
[~,Isort] = sort(ee,'descend');
I = I(Isort);
ee = round(ee(Isort));
for ii = 1:length(I)
  A = zeros(size(connMat{1,I(ii)}));
  nn = size(A,1);
  for jj = 1:6
    A = A+connMat{jj,I(ii)};
  end
  avgDens = nansum(A(:))/2 / (nn * (nn-1) /2) / jj;
  complexName = corum{I(ii),2};
  complexName = complexName(1:min([15 length(complexName)]));
  ss = sprintf('%s, %dx enriched, %1.2f avg density,  %d members,  %s',complexName,ee(ii),avgDens,corumN(I(ii)),corum{I(ii),4});
  disp(ss)
end



%% Make example connection matrices for enriched/not-enriched

I1 = 156; % 26S proteasome
I2 = 2693; % Emerin complex 25

figure
imagesc(connMat{5,I1})
colormap bone
xlabel('Protein number')
ylabel('Protein number')
axis square
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/connmat_enriched';
set(gcf,'paperunits','inches','paperposition',[.1 .1 4 4],'units','inches','position',[.1 .1 4 4])
print(sf,'-dpng');

figure
imagesc(connMat{4,I2})
colormap bone
xlabel('Protein number')
ylabel('Protein number')
axis square
sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/connmat_NOTenriched';
set(gcf,'paperunits','inches','paperposition',[.1 .1 4 4],'units','inches','position',[.1 .1 4 4])
print(sf,'-dpng');


%% Focus on a individual connection matricies

% show the connection matrices for ribosome from the different studies...
% We, and other CE-MS, find a lot of ribosome interactions.
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





%% Look at individual interactions

% full set of all interactions
A = [];
for ii = 1:length(interactomes)
  A = [A ;interactomes{ii}(:,3)];
end

NN = zeros(12,12);
JJ = zeros(12,12);
for ii = 1:length(interactomes)
  for jj = 1:length(interactomes)
    NN(ii,jj) = length(intersect(interactomes{ii}(:,3),interactomes{jj}(:,3)));
    UU = length(unique([interactomes{ii}(:,3) ;interactomes{jj}(:,3)]));
    JJ(ii,jj) =  NN(ii,jj) / UU;
  end
end



%% What's special about the 50?

% [~,II] = sort(nansum(inter2cor(1:3,:)),'descend');
% good1 = find(Ndatasets>2);
% II = II(ismember(II,good1));
% II = II(1:50);

II = I;


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
end

[unqID,I2] = unique(purID);
unqDesc = purDesc(I2);
pp = zeros(size(unqID));
EE = zeros(size(pp));
tmp2 = zeros(length(unqID),4);
for ii = 1:length(unqID)
  frac1 = [sum(purID(ID50==1)==unqID(ii)) length(purID(ID50==1))];
  frac2 = [sum(purID==unqID(ii)) length(purID)];
  tmp2(ii,:) = [frac1 frac2];
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



%% Write out the different corum subsets

% IMPORTANT!
% Keep apoptosis as validation set

fn = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Input/allComplexes.csv';
fid = fopen(fn);
head = fgetl(fid);

% S1: all studies
corumsub{1} = find(sum(EEp_complex([1:2 4:5 7:end],:) < pthresh) == 10);

% S2: all studies, but exlude the spliceosome
corumsub{2} = corumsub{1}(1:3);

% S3: almost all co-elution
corumsub{3} = find(sum(EEp_complex([1 2 4 5],:) < pthresh) >= 3);

% S4: all co-elution
corumsub{4} = find(sum(EEp_complex([1 2 4 5],:) < pthresh) == 4);

% S5: all co-elution, super stringent
corumsub{5} = find(sum(EEp_complex([1 2 4 5],:) < pthresh/100) == 4);

% S6: all co-elution, super-duper stringent
corumsub{6} = find(sum(EEp_complex([1 2 4 5],:) < pthresh/10^20) == 4);

% S7, S8: random
corumsub{7} = randsample(size(EEp_complex,2),length(corumsub{4}));
corumsub{8} = randsample(size(EEp_complex,2),length(corumsub{4}));

ww = corumN;
ww = ww .* (ww-1) / 2;

a = nanmean(ravg_corum(ds<7,:));
Rcor = nansum(a.*ww')/sum(ww);

disp(['average R, all of CORUM = ' num2str(Rcor)])
for ii = 1:length(corumsub)
  fn = ['/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/corumSubsets/allComplexes_sub' num2str(ii) '.csv'];
  fid = fopen(fn,'w');
  fprintf(fid,head);
  fprintf(fid,'\n');
  
  for jj = 1:length(corumsub{ii})
    I = corumsub{ii}(jj);
    for kk = 1:size(corum,2)
      fprintf(fid,'%s;',corum{I,kk});
    end
    fprintf(fid,'\n');
  end
  
  ravg = nansum(nanmean(ravg_corum(ds==3,corumsub{ii})) .* ww(corumsub{ii})') / sum(ww(corumsub{ii}));
  disp(['subset ' num2str(ii) ', average R = ' num2str(ravg)])
end

fclose all;



%% how many interactions per subset?

% craig, almost all, all, all stringent, all superstringent
%nintPerSub_50 = [6976 13677 22037 30248 70248 1398];
%nintPerSub_75 = [1080 2736 2936 3109 4723 0];
% craig, all, all stringent, all superstringent
nintPerSub_50 = [6976  22037 30248 70248 1398];
nintPerSub_75 = [1080  2936 3109 4723 0];

figure,hold on
bh = bar(nintPerSub_75);
set(bh,'FaceColor',[.1 .1 .1],'FaceAlpha',.2,'EdgeAlpha',1)
set(gca,'xtick',1:5,'xticklabel','','fontsize',11,'ytick',0:1000:5000)
ylabel('Number of predicted interactions','fontsize',12)

sf = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Justification document/Figure3.eps';
set(gcf,'paperunits','inches','paperposition',[.1 .1 5 2]*1.5,'units','inches','position',[.1 .1 5 2]*1.5)
print(gcf,sf,'-depsc2');



%% look at the corum subset used by Havugimana

fn = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/corumSubsets/Havugimana_subset.csv';
fid = fopen(fn);

fgetl(fid);
fgetl(fid);
fgetl(fid);

complexid = nan(10^4,1);
cc = 0;

while ~feof(fid)
  t1 = strsplit(fgetl(fid),',');
  ii = str2double(t1{1});
  ss = strrep(t1{2},'_',' ');
  ss = strrep(ss,'"', '');
  
  ss = strsplit(ss,'::');
  
  for jj = 1:length(ss)
    cc = cc+1;
    I = find(~cellfun('isempty',strfind(corum(:,2),ss{jj})));
    complexid(cc) = I(1);
  end
end
complexid = complexid(1:cc);
complexid = unique(complexid);



%% Is the overlap between replicates really 30%

clear fnn

fnn{1} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep1.csv';
fnn{2} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep2.csv';
fnn{3} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep3.csv';
fnn{4} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep4.csv';
fnn{5} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep5.csv';
fnn{6} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep6.csv';
fnn{7} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep7.csv';
fnn{8} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_Ap_rep8.csv';
% fnn{1} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep1.csv';
% fnn{2} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep2.csv';
% fnn{3} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep3.csv';
% fnn{4} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep4.csv';
% fnn{5} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep5.csv';
% fnn{6} = '/Users/Mercy/Academics/Foster/Jenny/mfiles/tmp/Nick_HeLa_rep6.csv';


ints = cell(1,length(fnn));
for ii = 1:length(fnn)
  tmp = readFinalInteractionList(fnn{ii});
  ints{ii} = tmp.text(:,1);
end

ints{1} = [ints{1}; ints{4}];
ints{2} = [ints{2}; ints{5}];
ints{3} = [ints{3}; ints{6}];
ints(4) = [];
ints(4) = [];
ints(4) = [];

JJ = nan(length(ints),length(ints));
for ii = 1:length(ints)
  for jj = 1:length(ints)
    JJ(ii,jj) = length(intersect(ints{ii},ints{jj})) / length(unique([ints{jj}; ints{ii}]));
  end
end


figure,imagesc(JJ),colorbar

