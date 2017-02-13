function check_prequisite(annot_fn, quic_dir)

    % check if dataset exists
    if(exist('dataset') ~= 2)
        throw(MException('twn:datasetnotexist', 'dataset() function does not exist in matlab.'))
    end

    % check if QUIC exists
    addpath(quic_dir);
    if(exist('QUIC') ~= 3)
        throw(MException('twn:quicnotexist', 'Either QUIC package (http://www.cs.utexas.edu/~sustik/QUIC/) was not installed properly or quic_directory setting was not set properly.'))
    end

    % check if QUIC runs properly
    try
        rng(101)
        M = rand(100,10);
        S = cov(M);
        [X W opt cputime iter_run dGap] = QUIC('default', S, 0.5, 1e-4, 1, 5);
    catch ME
        display(ME.message);
        throw(ME);
    end

    % check if transcript annotation file is OK
    try
        annot = dataset('File', annot_fn, 'ReadObsName', false, 'ReadVarNames', true, 'Delimiter', '\t');
        iso2gene = containers.Map(cellstr(annot(:,'transcript_id')), cellstr(annot(:,'gene_id')));
    catch ME
        display(ME.message);
        throw(ME);
    end
    
    disp('matlab installations are OK.');
end
