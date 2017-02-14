function twn(te_fn, iso_fn, out_fn, lambda_tt, lambda_ti, lambda_ii, lambda_d, lambda_s, iter, tol, verbose, standardize_data, annot_fn, quic_dir)

if ischar(lambda_tt)
    lambda_tt = str2num(lambda_tt);
end
if ischar(lambda_ti)
    lambda_ti = str2num(lambda_ti);
end
if ischar(lambda_ii)
    lambda_ii = str2num(lambda_ii);
end
if ischar(lambda_d)
    lambda_d = str2num(lambda_d);
end
if ischar(lambda_s)
    lambda_s = str2num(lambda_s);
end
if ischar(iter)
    iter = str2num(iter);
end
if ischar(tol)
    tol = str2num(tol);
end
if ischar(verbose)
    verbose = str2num(verbose);
end
if ischar(standardize_data)
    standardize_data = str2num(standardize_data);
end


addpath(quic_dir);

disp('reading inputs ...');
% skip header line, as gene names are modified when read in matlab
ds1 = dataset('File', te_fn, 'ReadObsNames', true, 'ReadVarNames', false, 'HeaderLines', 1, 'Delimiter', '\t'); 
samples1 = ds1.Properties.ObsNames;
M1 = double(ds1);

% read gene names
fh = fopen(te_fn);
first_line = fgets(fh);
genes = strsplit(first_line, {'\t', '\n'});
genes = genes(2:(size(genes,2)-1));
fclose(fh);


ds2 = dataset('File', iso_fn, 'ReadObsNames', true, 'ReadVarNames', false, 'HeaderLines', 1, 'Delimiter', '\t');
samples2 = ds2.Properties.ObsNames;
M2 = double(ds2);

% read isoform names
fh = fopen(iso_fn);
first_line = fgets(fh);
isoforms = strsplit(first_line, {'\t', '\n'});
isoforms = isoforms(2:(size(isoforms,2)-1));
fclose(fh);

features = [genes, isoforms];
[nse ne] = size(M1);
[nsi ni] = size(M2); 
nfeatures = ne + ni;

%% error if the samples are not in the same order
if nse ~= nsi
    throw(MException('twn:diffnoofsample', 'total expression and isoform ratio data must have same no. of samples.'))
end

for i = 1:size(samples1,1)
    if ~strcmp(samples1{i}, samples2{i})
        throw(MException('twn:sampleorder', 'samples must maintain the same order in both expression and isoform data'))
    end
end


M = [M1(:,1:ne), M2(:,1:ni)];

% standardize data
if standardize_data==1 
    M = standardize(M')';
end


%% read transcript annotation file and process it.
annot = dataset('File', annot_fn, 'ReadObsName', false, 'ReadVarNames', true, 'Delimiter', '\t');
iso2gene = containers.Map(cellstr(annot(:,'transcript_id')), cellstr(annot(:,'gene_id')));
gene2isoidx = containers.Map(); 
for idx = 1:size(isoforms,2)
    iso = char(isoforms(idx));
    gene = char(iso2gene(iso));
    try
        gene2isoidx(gene) = [gene2isoidx(gene) idx];
    catch
        gene2isoidx(gene) = [idx];
    end
end


disp('creating penalty metrix ...');
% create lambda matrix
L = zeros(nfeatures, nfeatures);
L(1:ne, 1:ne) = lambda_tt;
L(1:ne, ne+1:nfeatures) = lambda_ti;
L(ne+1:nfeatures, 1:ne) = lambda_ti;
L(ne+1:nfeatures, ne+1:nfeatures) = lambda_ii;
L(logical(eye(nfeatures))) = lambda_d; % diagonal penalty
% put small lambda between isoforms of same geme
for k = gene2isoidx.keys()
    indexes = gene2isoidx(char(k));
    if size(indexes,2) <= 1
        continue;
    end
    
    for i = 1:(size(indexes,2)-1)
        idx1 = indexes(i);
        for j=(i+1):size(indexes,2)
            idx2 = indexes(j);
            L(ne+idx1, ne+idx2) = lambda_s;
            L(ne+idx2, ne+idx1) = lambda_s;
        end
    end
end
% put small lambda between te and isoforms of same geme
for idx1 = 1:size(genes,2)
    g = genes(1,idx1);
    try
        indexes = gene2isoidx(char(g));
    catch
        continue;
    end

    if size(indexes,2) < 1
        continue;
    end

    for j=1:size(indexes,2)
        idx2 = indexes(j);
        L(idx1, ne+idx2) = lambda_s;
        L(ne+idx2, idx1) = lambda_s;
    end
end

% run quic
disp('running QUIC');
S = cov(M);
[X W opt cputime iter_run dGap] = QUIC('default', S, L, tol, verbose, iter);

clear S;
clear L;


% save sparse X
disp('saving QUIC results');
display(datestr(now));
[rows,cols,values] = find(X);
clear X;
to_take = rows<=cols;
rows_taken = rows(to_take);
cols_taken = cols(to_take);
name1s = features(rows_taken);
name2s = features(cols_taken);
weights = values(to_take);
types = ones(size(rows_taken,1),1);
types(rows_taken <= ne & cols_taken > ne) = 2;
types(rows_taken > ne & cols_taken > ne) = 3;

quic_out_fn = sprintf('%s.quic.txt', out_fn);
fid = fopen(quic_out_fn, 'w');
fprintf(fid, 'Name1\tName2\tEdge type\tEdge weight\n');
%strings= strcat(name1s,{'\t'}, name2s, {'\t'}, strread(num2str(types'),'%s')', {'\t'},strread(num2str(weights'),'%s')', {'\n'}));
C = reshape( cat(4, name1s, name2s, strread(num2str(types'),'%s')', strread(num2str(weights'),'%s')'), size(types,1),4);
formatSpec = '%s\t%s\t%s\t%s\n';
for row = 1:size(types,1)
    fprintf(fid,formatSpec,C{row,:});
end
fclose(fid);
display(datestr(now));


% save QUIC outputs
quic_info_out_fn = sprintf('%s.quic.info', out_fn);
fid = fopen(quic_info_out_fn, 'w');
fprintf(fid, 'optimum objective: %f\n', opt);
fprintf(fid, 'iteration: %d\n', iter_run);
fprintf(fid, 'time to run: %f\n', cputime);
fclose(fid);


