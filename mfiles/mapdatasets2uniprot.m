
clear fn mf
fn{1} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/Bandyopadhyay2010_pairs.csv';
mf{1} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/Bandyopadhyay2010_map.fasta';

fn{2} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/ewing2007_pairs.csv';
mf{2} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/ewing2007_map.fasta';

fn{3} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/rual2005_pairs.csv';
mf{3} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/rual2005_map.fasta';

fn{4} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/stelzl2005_pairs.csv';
mf{4} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/stelzl2005_map.fasta';

fn{5} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/wang2011_pairs.csv';
mf{5} = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/wang2011_map.fasta';

fnout = {'/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Bandyopadhyay_et_al.csv'...
  '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Ewing_et_al.csv'...
  '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Rual_et_al.csv'...
  '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Stelzl_et_al.csv'...
  '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot/Wang_et_al.csv'};


for ii = 1:5
  data = fastaread(mf{ii});
  gene2prot = cell(length(data),2);
  for jj = 1:length(data)
    t1 = strsplit(data(jj).Header,'|');
    gene2prot{jj,2} = t1{2};
    
    tmp = strsplit(t1{3},'GN=');
    tmp2 = strsplit(tmp{2},' ');
    gene2prot{jj,1} = tmp2{1};
  end
  
  fid = fopen(fn{ii});
  fidout = fopen(fnout{ii},'w');
  while ~feof(fid)
    t1 = strsplit(fgetl(fid),',');
    [ia1,ib] = find(ismember(gene2prot,t1{1}));
    [ia2,ib] = find(ismember(gene2prot,t1{2}));
    if ~isempty(ia1) && ~isempty(ia2)
      fprintf(fidout,'%s,%s\n',gene2prot{ia1,2},gene2prot{ia2,2});
    end
  end
  fclose all;
end


%% convert stelzl_et_al.rtf to stelzl_et_al.csv

fn = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/stelzl_et_al.rtf';
fid = fopen(fn);
cc = 1;
stelzl = cell(10^5,2);
while ~feof(fid)
  t1 = strsplit(fgetl(fid),' ');
  cc1 = str2double(t1{1});
  if cc==cc1
    stelzl{cc,1} = t1{3};
    stelzl{cc,2} = t1{6};
    cc = cc+1;
  end
end
stelzl = stelzl(1:cc1,:);
fclose(fid);

fn = '/Users/Mercy/Academics/Foster/PCP-SILAC/No-reference methods/Systematic bias?/uniprot_badformat/stelzl2005_pairs.csv';
fid = fopen(fn,'w');
for ii = 1:cc1
  fprintf(fid,'%s,%s\n',stelzl{ii,1},stelzl{ii,2});
end
