
% ************************************************************************
%
%   PLISH Probe Design Program, Version "BLAST210", 02/2020
%    
%   By. A. M. Andruska
%
%   Please credit the author if results end in publication.
%
%*************************************************************************
clear all;
close all;

% UPDATES FOR "BLAST210" VERSION (02/2020):
% Now program uses BLAST 2.10+ which has new structure for filtering by
% taxon. Instead of .gi files, user inputs TaxID in a dialog box when
% running the script. 
% Example Taxids:Variants 
%   Homo Sapiens:               9606
%   Mus musculus:               10090
%   Rattus norvegicus:          10116
%   Caenorhabditis elegans:     6239     
%
% This version uses the V5 version of the NCBI refseq_rna database

% UPDATES FOR "Santa Cruz" VERSION (07/2018):
% This is a combination of programs D2 and R2. It will create an easier to
% use output program that bases successful target sites on the score. The
% socre cutoff is adjustable, but a cutoff of 30 is what matches the cutoff
% on the online "Human G+T" database and seems reasonable.


% How to use:
%   1. Pick an mRNA and copy the entire sequence (numbers and all) from
%   NCBI nucleotide site.
%   2. This program requires that two species files are present in the same
%   directory as the program:
%       1. mus.gi.txt
%       2. homosapiens.gi.txt
%   3. Run the program in MATLAB
%   4. The gene name given is a unique label that helps to identify files
%   that the program uses, in addition to labeling the H-probes.
%   5. Open the output file in opera (works better than notepad)
%       Output forms:
%           <GeneName>_<Fluorophore>_Output_File_<Date>_All_Blast_Results.txt
%           <GeneName>_<Fluorophore>_Program_Output_<Date>_With_Hits_Limited_To_<Hit Limit>.txt

% MODIFICATIONS THAT NEED TO BE MADE TO MATLAB FUNCTIONS IN ORDER TO MAKE
% THE PROGRAM WORK CORRECTLY:
%
% Modify the function "blastreadlocal". 
% On line 492, due to the way the blastoutput file is wrtten, I changed the MATLAB code from
%
%   if strfind(blasttext,'Strand =') % Strand is included with nucleotide sequences
%       to
%   if strfind(blasttext,'Strand=') % Strand is included with nucleotide sequences


% USER SPECIFIED OPTIONS HEADER:
% ------------------------------------------------------------------------

% 1 "mode_option" The specific site at which to design a probe. Options:
    %       1 = 5' - A|CC - 3' with exclusion of GGT/GGC/GGG/GAT before*
    %       2 = 5' - A|CC - 3' with no preceeding exclusion criteria
    %       3 = 5' - A|G - 3' or 5' - T|A - 3' Version
    % 
    % * This is the presumed optimal Holiday Junction Configuration
    %
    % This is reflected in the output file as "_Accplus.txt" for 1,
    % "_Accminus.txt" for 2, and "_TAAG.txt" for 3
    % mode_option = 3;
    % This is now a user specified dialog box
    
% 2 "fluoro_option" Will specifiy the tail to which the probes are congugated
    %       Left:   5' - LLLLL:hprobe - 3'
    %       Right:  5' - hprobeRC:RRRRR - 3'
    %
    %   Current Options:
    %       1 = Cy3 / HLR3X
    %       2 = Cy5 / HLR5X
    %       3 = a488 / HLR2X
    %       4 = Texas Red / HLR4X
    %       5 = PB405 / HLR6X
    %
    % NOTE: This has now been replaced by a seperate menu which allows one
    % to choose from all of the current bridges at the start of the
    % program.
    
% 3 "species_option" 
    %       1 = human
    %       2 = mouse
    % This is now directly set with a message box that pops up.
    
% 4 "hlength" specifies the total BPs the H-probes will cover (we use 40)
    hlength = 40;
    
% 5 "dimer_length" sets the number of base pairs which must align in order
    % to be considered a dimer
    dimer_length = 6;
    
% 6 "library_directory"
% This is a directory to which each new probe is saved. 
    library_directory = 'C:\Users\yalbi\OneDrive\Documents\GitHub\PLISH_ProbeDesigner\hprobelibrary';
    
% 7 "Hit Score Cutoff"
% This is the criteria by which we say that an alignment is significnat or
% not. This has been set to 30, but could be adjusted to a different score.
% This is now an input dialog box, but this serves as the default.
    bitscore_sig = 30;
    
% 8 "List Length"
% This is how many probes will be put in the final output list.
    numtoprint = 250;
    
% 9 "Enviromental Variables"
% Need to specify the enviromental variables that the local NCBI blast
% program uses:
%setenv('PATH', 'C:\NCBI\blast-2.10.0+\bin');
%setenv('BLASTDB', 'C:\NCBI\blast-2.10.0+\db'); 
    
% NCBI BLAST USER SPECIFIED OPTIONS:
%--------------------------------------------------------------------------

% 1 Database:   X
    % refseq_rna contains all RNA sequences which start as NM_, NR_, XM_,
    % and XR_ (there will be multiple species returned)
    % This is installed under C:\NCBI\blast-2.10.0+\db
    database_ncbi = 'refseq_rna';
    
% 2 Entrez Query:
    % This string will help limit to species of interest... it is
    % automatically set below by the "species_option" flag
    % This is not supported on local searches, only remote searches.
    
% 3 Maximum number of hits to return. Default is 250, now set to 100.
    maxseq_ncbi = 100;
    
% 4 Expect:     X
    % Statistical significance threshold for matches against database sequences
    % Default is 10, we have set 100 for more stringency
    expect_ncbi = 100;
    
% 5 Word:       X
    % Word length for the query sequence. Options are 7, 11, 15
    word_ncbi = 7;
    
% 6 Match scores. Matching and mismatching scores in a nucleotide alignment, 
    % This sets the Reward for a nucleotide match = "ncbi_match"
    % This sets the Reward for a nucleotide mismatch  = "ncbi_mismatch"
    % Default is 1 and -2 for blastn
    % Online, these are set to 2 and -3
    ncbi_match = 2;
    ncbi_mismatch = -3;
    
% 7 Gap Scores; cost to create and extend a gap in an alignment. 
    % Default on blastn is 5 and 2, so this is not changed here
    % ncbi_gap = [5,2];
    
% 8 Program. This is generally always "blastn" but can be changed to others.
    program_ncbi = 'blastn';
    
    
% The below options are attempts at pairing down the probes to successful
% sequences:

     
% 09 NEW: Score cutoff:
% If a hit has a score less than this then it is considered to be OK.
     bitscore_sig = 30;

% 10 BITSCORE Range:
% The probes are cycled over this range of bitscores to get an idea of how
% good they are
    bitscore_low = 10;
    bitscore_high = 30;

% END OF USER SPECIFIED OPTIONS.
%--------------------------------------------------------------------------


% Part 00: Start from Scratch, pre-load data, or make a printout of probes
load_option = 1;
option_list = {    'Start from Scratch with a new Gene',...                                     % 1
            'Load an Existing Gene that has already been run',...                               % 2
            'Make an Excel Sheet of Probes for a gene that has already been run'};              % 3
[load_option,tf] = listdlg(    'ListString',option_list,...
                                'SelectionMode','single',...
                                'ListSize',[450,250],...
                                'Name','Select the Mode of Operation',...
                                'OKString','Apply',...
                                'CancelString','Default (Start from Scratch)');


if load_option == 1


    % Part 0: Ask for an input for all the transcript Variants:
    tv_cell = {};
    while isempty(tv_cell)
        tvs = inputdlg('Please enter the Accession numbers of all transcript Variants separated by a Return:',...
            'Enter List of Transcript Variant Accession Numbers (IE  NM_001204.6)', [10 40]);
        tvs = tvs{1,1};
        no_tvs = size(tvs,1);
        for i = 1:no_tvs
            tv_cell = [tv_cell; lower(strtrim(tvs(i,:))) ];
        end
    end

    % Part 1: Input the mRNA sequence and clean it up
    tmRNA = inputdlg('Please paste the entire sequence with index numbers below:',...
        'Enter mRNA Gene Sequence from NHBI', [40 100]);

    gene_name = inputdlg('Enter The Gene Name',...
                 'Name Entry', [1 50]);
    gene_name = gene_name{1,1};

    temp_bit_score = inputdlg('Choose the BitScore ("Score") that defines appropriate homology. A score of 30 or less is the default.',...
                 'Bitscore', [1 50]);
    temp_bit_score = str2num(temp_bit_score{1,1});
    if temp_bit_score ~= bitscore_sig
        bitscore_sig = temp_bit_score;
    end
    
    tmRNA = tmRNA{1,1};
    s_tmRNA = size(tmRNA);

    % This line will make a new directory to store this data
    savepath = [library_directory,'\',gene_name,'_ver_02_2020'];
    status = mkdir(savepath);

    % Part 1-B: Ask for an input regarding the taxid limitation
    taxid_temp = {};
    ncbi_taxid = '9606';
    while isempty(taxid_temp)
        taxid_temp = inputdlg('Please enter the NCBI Taxid to limit the BLAST Database. Homo sapiens is 9606. Mus musculus is 10090.:',...
            'Enter Taxid ID Number', [10 40]);
        ncbi_taxid = (taxid_temp{1,1});
    end
    
    mRNA = '';
    for i = 1:s_tmRNA(1,1)
        for j = 1:s_tmRNA(1,2)
            q = str2num(tmRNA(i,j));
            if isempty(q)
                mRNA = strcat(mRNA,tmRNA(i,j));
            end
        end
    end  



    % Part 2: All the variable bridges for future use are coded in here. This is somewhat
    % cumbersome, but will be stored in a structue which will make it easier to
    % use later. 
    vbseqs = struct();
    for i = 1:28
        if i == 1
                % VB1
                comp_l = 'AGGTCAGGAATACTTACGTCGTTATGG';       % HL3X:LLLL
                comp_r = 'TTATAGGTCGAGTAGTATAGCCAGGTT';       % RRRR:HR3X
                pnl = 'HLC2-VB01-';
                pnr = 'HRC2-VB01-';     
                id = 'VB1';
        elseif i == 2
                % VB2
                comp_l = 'AGGTCAGGAATACTTAGCTATTGATGG';
                comp_r = 'TTCTGTGTAGACGACTATAGCCAGGTT';
                pnl = 'HLC2-VB02-';
                pnr = 'HRC2-VB02-'; 
                id = 'VB2';
        elseif i == 3
                % VB3
                comp_l = 'AGGTCAGGAATACCAGGTTGTAATGG';
                comp_r = 'TTAGTATGATGACAATATAGCCAGGTT';
                pnl = 'HLC2-VB03-';
                pnr = 'HRC2-VB03-'; 
                id = 'VB3';
        elseif i == 4
                % VB4
                comp_l = 'AGGTCAGGAATACGACTACGAGTAGG';
                comp_r = 'TTCTTGTGGTGTAAGTATAGCCAGGTT';
                pnl = 'HLC2-VB04-';
                pnr = 'HRC2-VB04-'; 
                id = 'VB4';
        elseif i == 5
                % VB5
                comp_l = 'AGGTCAGGAATACGTGGAGTGACCCGG';
                comp_r = 'TTACCCCTTGTGAGATATAGCCAGGTT';
                pnl = 'HLC2-VB05-';
                pnr = 'HRC2-VB05-'; 
                id = 'VB5';
        elseif i == 6
                % VB6
                comp_l = 'AGGTCAGGAATAAGACATGCTGACAGG';
                comp_r = 'TTCACACAATATAGATATAGCCAGGTT';
                pnl = 'HLC2-VB06-';
                pnr = 'HRC2-VB06-'; 
                id = 'VB6';
        elseif i == 7
                % VB7
                comp_l = 'AGGTCAGGAATATTCATATTGGTACGG';
                comp_r = 'TCCAGCATGGAAGATTATAGCCAGGTT';
                pnl = 'HLC2-VB07-';
                pnr = 'HRC2-VB07-'; 
                id = 'VB7';
        elseif i == 8
                % VB8
                comp_l = 'AGGTCAGGAATAGGACAACAAGGTCGG';
                comp_r = 'TCAGCGCTAATCACATATAGCCAGGTT';
                pnl = 'HLC2-VB08-';
                pnr = 'HRC2-VB08-';
                id = 'VB8';
        elseif i == 9
                % VB9
                comp_l = 'AGGTCAGGAATACACTGGGCACGGAGG';
                comp_r = 'TCATGAAGAAAAGAATATAGCCAGGTT';
                pnl = 'HLC2-VB09-';
                pnr = 'HRC2-VB09-'; 
                id = 'VB9';
        elseif i == 10
                % VB10
                comp_l = 'AGGTCAGGAATAGTAAACAACCCATGG';
                comp_r = 'TGCGAGAGCCGGAATTATAGCCAGGTT';
                pnl = 'HLC2-VB10-';
                pnr = 'HRC2-VB10-'; 
                id = 'VB10';
        elseif i == 11
                % VB11
                comp_l = 'AGGTCAGGAATAATCGGAAGGCAGTGG';
                comp_r = 'TTGAGGCTGCTTATCTATAGCCAGGTT';
                pnl = 'HLC2-VB11-';
                pnr = 'HRC2-VB11-'; 
                id = 'VB11';
        elseif i == 12
                % VB12
                comp_l = 'AGGTCAGGAATAGTCTCAGCTCCGAGG';
                comp_r = 'TGGATTTACTTCCAGTATAGCCAGGTT';
                pnl = 'HLC2-VB12-';
                pnr = 'HRC2-VB12-'; 
                id = 'VB12';
        elseif i == 13
                % VB13
                comp_l = 'AGGTCAGGAATATAGGCCGAAAAAAGG';
                comp_r = 'TTGTATACTATTAAATATAGCCAGGTT';
                pnl = 'HLC2-VB13-';
                pnr = 'HRC2-VB13-'; 
                id = 'VB13';
        elseif i == 14
                % VB14
                comp_l = 'AGGTCAGGAATAGTGGCTTGTTATGGG';
                comp_r = 'TGAATGGACTAGCACTATAGCCAGGTT';
                pnl = 'HLC2-VB14-';
                pnr = 'HRC2-VB14-'; 
                id = 'VB14';
        elseif i == 15
                % VB15
                comp_l = 'AGGTCAGGAATACGGTGTTCGTTTTGG';
                comp_r = 'TCCGCTGTTCTACACTATAGCCAGGTT';
                pnl = 'HLC2-VB15-';
                pnr = 'HRC2-VB15-'; 
                id = 'VB15';
        elseif i == 16
                % VB16
                comp_l = 'AGGTCAGGAATAGACTACTCAAACTGG';
                comp_r = 'TGATTACACCTTCAGTATAGCCAGGTT';
                pnl = 'HLC2-VB16-';
                pnr = 'HRC2-VB16-'; 
                id = 'VB16';
        elseif i == 17
                % VB17
                comp_l = 'AGGTCAGGAATATGTCAATTGATGGGG';
                comp_r = 'TTAAGAGAACTTTCTTATAGCCAGGTT';
                pnl = 'HLC2-VB17-';
                pnr = 'HRC2-VB17-'; 
                id = 'VB17';
        elseif i == 18
                % VB18
                comp_l = 'AGGTCAGGAATAAATTTTAGCCAATGG';
                comp_r = 'TAGCATATCTTCTATTATAGCCAGGTT';
                pnl = 'HLC2-VB18-';
                pnr = 'HRC2-VB18-'; 
                id = 'VB18';
        elseif i == 19
                % VB19
                comp_l = 'AGGTCAGGAATATGACTTTTCCGCGGG';
                comp_r = 'TCTTTTTGTAACAATTATAGCCAGGTT';
                pnl = 'HLC2-VB19-';
                pnr = 'HRC2-VB19-'; 
                id = 'VB19';
        elseif i == 20
                % VB20
                comp_l = 'AGGTCAGGAATACTTCGCGAAACCAGG';
                comp_r = 'TATAGTAACGAGGACTATAGCCAGGTT';
                pnl = 'HLC2-VB20-';
                pnr = 'HRC2-VB20-'; 
                id = 'VB20';
        elseif i == 21
                % VB21
                comp_l = 'AGGTCAGGAATAAAGGGCCCCTGCCGG';
                comp_r = 'TGCATCCGCACGATTTATAGCCAGGTT';
                pnl = 'HLC2-VB21-';
                pnr = 'HRC2-VB21-'; 
                id = 'VB21';
        elseif i == 22
                % VB22
                comp_l = 'AGGTCAGGAATAAGCGACAAGTTTGGG';
                comp_r = 'TAGGAGCGACGACGATATAGCCAGGTT';
                pnl = 'HLC2-VB22-';
                pnr = 'HRC2-VB22-'; 
                id = 'VB22';
        elseif i == 23
                % VB23
                comp_l = 'AGGTCAGGAATATGTCAGGCGAGGCGG';
                comp_r = 'TTTATACACGAGACTTATAGCCAGGTT';
                pnl = 'HLC2-VB23-';
                pnr = 'HRC2-VB23-'; 
                id = 'VB23';
        elseif i == 24
                % VB24
                comp_l = 'AGGTCAGGAATACTGTACATACGGGGG';
                comp_r = 'TCCTAACCGTGCGCCTATAGCCAGGTT';
                pnl = 'HLC2-VB24-';
                pnr = 'HRC2-VB24-'; 
                id = 'VB24';
        elseif i == 25
                % Cy3
                comp_l = 'TATTCGTTCGAACTTACGTCGTTATG';       % HL3X:LLLL
                comp_r = 'TTATACGTCGAGTTGACCGACGTATTG';      % RRRR:HR3X
                pnl = 'HL3X-';
                pnr = 'HR3X-';   
                id = 'Cy3';
        elseif i == 26
                comp_l = 'TCGTACGTCTAACTTACGTCGTTATG';      % HL2X
                comp_r = 'TTATACGTCGAGTTGAAGAACAACCTG';     % HR2X
                pnl = 'HL2X-';
                pnr = 'HR2X-';
                id = 'a488';
        elseif i == 27
                comp_l = 'TTAGTAGGCGAACTTACGTCGTTATG';    % HL4X
                comp_r = 'TTATACGTCGAGTTGAACATAAGTGCG';   % HR4X  
                pnl = 'HL4X-';
                pnr = 'HR4X-';
                id = 'TxRd';
        else
                comp_l = 'TAGGTCAGGAAACTTACGTCGTTATG';      % HL6X
                comp_r = 'TTATACGTCGAGTTGAATAGCCAGGTT';     % HR6X
                pnl = 'HL6X-';
                pnr = 'HR6X-';
                id = 'Cy5';
        end

        vbseqs(i).comp_l = comp_l;
        vbseqs(i).comp_r = comp_r;
        vbseqs(i).pnl = pnl;
        vbseqs(i).pnr = pnr;
        vbseqs(i).name = id;
    end


    h = waitbar(0,'1','Name','Designing H-Probes...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0)

    % Part 3 Big loops to find and create the hprobe data structure
    h_probe_data = struct();
    total_length = (length(mRNA)-21)-21;
    fastafilestring = [gene_name, '_fasta_file.txt'];
    counter = 1;
    waitbarcount = 3*total_length;
    waitbarcount2 = 0;

    for m = 1:3

        % Conditions from which to determine the probe location
        if m == 1
            tseq1 = 'acc';
            tseq2 = 'hhh';
            check_for_bad_preceeding = 1;
            tseq1_l = 2;
            tseq2_l = 2;
            excel_name = '_Accplus';
        elseif m == 2
            tseq1 = 'acc';
            tseq2 = 'hhh';
            check_for_bad_preceeding = 0;
            tseq1_l = 2;
            tseq2_l = 2;
            excel_name = '_Accminus';
        elseif m == 3
            tseq1 = 'ta';
            tseq2 = 'ag';
            check_for_bad_preceeding = 0;
            tseq1_l = 1;
            tseq2_l = 1;
            excel_name = '_TAAG';
        end

        for i = 21:(length(mRNA)-21)

            waitbarcount2 = waitbarcount2 + 1;

            waitbar( (i/(length(mRNA)-21)), h, strcat( num2str( (waitbarcount2/waitbarcount)*100, '%12.1f'),['% Complete, Step ',num2str(m),' of 3']) );

            if strcmp(mRNA(1,i:(i+tseq1_l)),tseq1) || strcmp(mRNA(1,i:(i+tseq2_l)),tseq2)
                % Identify the bp# that identifies this potential oligo
                oligo_bp_start = (i-19);
                oligo_l = mRNA(1, (i-19):i);
                oligo_r = mRNA(1, (i+1):(i+20));

                % This figures out bp 17-19, which should not be {GGT, GGC,
                % GAT, or GGG}
                bad_seq = oligo_l(1, 17:19);
                bad_flag = 0;

                if check_for_bad_preceeding == 1
                    if strcmp('GGT', bad_seq)
                        bad_flag = 1;
                    elseif strcmp('GGC', bad_seq)
                        bad_flag = 1;
                    elseif strcmp('GAT', bad_seq)
                        bad_flag = 1;
                    elseif strcmp('GGG', bad_seq)
                        bad_flag = 1;
                    end
                end

                if bad_flag == 0
                    % At the point we've found an acceptable h-probe...
                    h_probe_data(counter).gene = gene_name;
                    h_probe_data(counter).bp_start = i-19;
                    h_probe_data(counter).design_mode = m;
                    h_probe_data(counter).mRNA_sequence = [oligo_l, oligo_r];
                    h_probe_data(counter).HL_Seq_20 = seqrcomplement(oligo_l);
                    h_probe_data(counter).HR_Seq_20 = seqrcomplement(oligo_r);
                    h_probe_data(counter).threeprimecloseness = 100*( (i-19)/length(mRNA) );

                    warnState = warning; %Save the current warning state
                    warning('off','Bioinfo:fastawrite:AppendToFile');
                    fastawrite(fastafilestring, num2str(counter), [oligo_l, oligo_r] );
                    warning(warnState); %Reset warning state to previous settings    

                    blastoutputfilestring = [gene_name,'_blastoutputfile.txt'];
                    dosCommand = [  'blastn -query ', fastafilestring, ...
                                    ' -taxids ', ncbi_taxid, ...
                                    ' -db ',database_ncbi, ...
                                    ' -evalue ', num2str(expect_ncbi),...
                                    ' -word_size ', num2str(word_ncbi),...
                                    ' -gapopen 5',...
                                    ' -gapextend 2',...
                                    ' -dust no',...
                                    ' -num_alignments ', num2str(maxseq_ncbi),...
                                    ' -reward ', num2str(ncbi_match), ...
                                    ' -penalty ', num2str(ncbi_mismatch), ...
                                    ' -out ',blastoutputfilestring];
                    [status,cmdout] = dos(dosCommand);
                    h_probe_data(counter).ncbi = blastreadlocal(blastoutputfilestring,0);
                    delete(fastafilestring);
                    delete(blastoutputfilestring);

                    ncbi_short_temp = [];

                    number_tvs_included = 0;
                    % This section will parse the blast data output and
                    % construct the score matrices which are used for
                    % generating reports later on
                    for j = 1:size(h_probe_data(counter).ncbi.Hits,2)
                        name = h_probe_data(counter).ncbi.Hits(j).Name;
                        strand = h_probe_data(counter).ncbi.Hits(j).HSPs.Strand;
                        score = h_probe_data(counter).ncbi.Hits(j).HSPs.Score;
                        expect = h_probe_data(counter).ncbi.Hits(j).HSPs.Score;

                        tvflag = 0;
                        plus_plus_flag = 0;

                        % Identify if the hit is important
                        if strcmp(strand,'Plus/Plus')

                            % Test if the hit is a transcript Variant
                            for k = 1:size(tv_cell, 1)
                                if contains(lower(name), tv_cell{k,1})
                                    tvflag = 1;
                                end
                            end
                            plus_plus_flag = 1;
                        end

                        ncbi_short_temp = [ ncbi_short_temp, [ j; plus_plus_flag; score; tvflag ] ];
                    end

                    h_probe_data(counter).ncbi_short = ncbi_short_temp;   
                    h_probe_data(counter).ncbi_short_legend = {'Alignment Number'; 'Plus/Plus (1) or Not (0)'; 'Bitscore';'Transcript Variant (1) or not (0)'};


                    % Narrow to only Plus/Plus Results
                    if isempty(ncbi_short_temp) == 0
                        ncbi_report = ncbi_short_temp(:, find(ncbi_short_temp(2,:) == 1));
                    else
                        off_target_hits = 0;
                    end
                    if isempty(ncbi_report) == 0
                        % Narrow further to only the bitscore
                        ncbi_report = ncbi_report(:, find(ncbi_report(3,:) >= bitscore_sig));    
                        % Exclude the transcript Variants to get the total number
                        % of significant off target hits for the given probe
                        off_target_hits = size( ncbi_report(:, find( ncbi_report(4,:) == 0 )), 2);
                    else
                        off_target_hits = 0;
                    end

                    % Update the blast reports
                    h_probe_data(counter).ncbi_report = ncbi_report;
                    h_probe_data(counter).total_off_target_hits = off_target_hits;

                    % Include data for transcript Variants
                    if isempty(ncbi_short_temp) == 0
                        h_probe_data(counter).tvs_targeted = sum(ncbi_short_temp(4,:));
                    else
                        h_probe_data(counter).tvs_targeted = 0;
                    end
                    h_probe_data(counter).total_tvs = tv_cell;
                    h_probe_data(counter).significant_bitscore = bitscore_sig;


                    % This section will now figure out all of the thermodynamic
                    % parameters of the oligo after it has been conjugated with
                    % the proper variable bridge

                    vb = struct();

                    for j = 1:size(vbseqs,2)

                        temp_left = [ vbseqs(j).comp_l, seqrcomplement(oligo_l) ];
                        temp_right = [ seqrcomplement(oligo_r), vbseqs(j).comp_r ];          
                        temp_left_name = [ vbseqs(j).pnl, num2str(m), '-', gene_name, '-', num2str(i-19) ];
                        temp_right_name = [ vbseqs(j).pnr, num2str(m), '-', gene_name, '-', num2str(i-19) ];

                        h_p_l_prop = oligoprop(temp_left, 'Dimerlength', dimer_length);
                        h_p_r_prop = oligoprop(temp_right, 'Dimerlength', dimer_length);

                        % Build the results into 
                        vb(j).Full_Left_Probe = temp_left;
                        vb(j).Left_Name = temp_left_name;
                        vb(j).Full_Right_Probe = temp_right;
                        vb(j).Right_Name = temp_right_name;
                        vb(j).Thermo_Left = h_p_l_prop;
                        vb(j).Thermo_Right = h_p_r_prop;
                        vb(j).Tm_Left = mean(h_p_l_prop.Tm);
                        vb(j).Tm_Right = mean(h_p_r_prop.Tm);
                        vb(j).No_Dimers_Left = size(h_p_l_prop.Dimers,1);
                        vb(j).No_Dimers_Right = size(h_p_r_prop.Dimers,1);
                        vb(j).No_Hairpins_Left = size(h_p_l_prop.Hairpins,1);
                        vb(j).No_Hairpins_Right = size(h_p_r_prop.Hairpins,1);
                        vb(j).Delta_G_Left = mean(h_p_l_prop.Thermo(:,3));
                        vb(j).Delta_G_Right = mean(h_p_r_prop.Thermo(:,3));
                        vb(j).MolWt_Left = h_p_l_prop.MolWeight;
                        vb(j).MolWt_Right = h_p_r_prop.MolWeight;

                    end
                    h_probe_data(counter).final_h_probes = vb;

                    counter = counter + 1;
                end
            end
        end
    end

    delete(h);

    % This section comes out with a matrix that can be used to sort the probes
    % base on blast results and transcript Variant hits

    sorting_matrix = [];
    for i = 1:size(h_probe_data,2)
        sorting_matrix = [sorting_matrix; [i, h_probe_data(i).total_off_target_hits, h_probe_data(i).tvs_targeted, h_probe_data(i).bp_start ] ];
    end
    ranking_matrix = sortrows(sorting_matrix, [2,3],{'ascend','descend'} );


 
    % This section will vary the bitscore from a low number to a high
    % number. It will then count up the nubmer of times each probe passes
    % (IE has no significant hits) for a given bitscore. It first has to
    % get rid of any duplicate probesites
    probesites = [-1];
    probeindicies = [-1];
    for i = 1:size(h_probe_data,2)
        tbp = h_probe_data(i).bp_start;
        if ismember(probesites, tbp) == 0
            probesites = [probesites, tbp];
            probeindicies = [probeindicies, i];
        end
    end
    probeindicies = probeindicies(1, 2:size(probeindicies,2));
    probesites = probesites(1, 2:size(probesites,2));
    % Now cycle over bitscore and probesite to find out if the probes are
    % good or not
    probe_hits_sum = zeros(size(probesites));
    for i = bitscore_low:bitscore_high
        for j = 1:size(probesites,2)
            tn = 0;
            temp_mat = h_probe_data(j).ncbi_short;
            if isempty(temp_mat) == 0
                temp_mat = temp_mat(:, find(temp_mat(4,:) == 0) );
                if isempty(temp_mat) == 0
                    temp_mat = temp_mat(3, find(temp_mat(2,:) == 1) );
                    tn = length(find(temp_mat > i));
                end
            end
            if tn > 0
                probe_hits_sum(1,j) =  probe_hits_sum(1,j) + 1;
            end
        end
    end
    resort_mat = [probesites',probeindicies',probe_hits_sum'];
    resort_mat = sortrows(resort_mat,[3],'ascend');

    
    % At this point save the data structure and the mRNA
    filenamemat = [savepath,'\',gene_name,'_DataStructure_',date,'.mat'];
    save(filenamemat, 'h_probe_data','mRNA','ranking_matrix','tv_cell','resort_mat','vbseqs','dimer_length');
    
    
   
    
    
    
    % This will make a figure showing the probe sites, and respective homology
    plot_matrix = sortrows(ranking_matrix, [4],{'ascend'} );
    max_off_target = max( plot_matrix(:,2) );
    max_tvs = size(h_probe_data(1).total_tvs, 1);
    off_t_plot = plot_matrix(:,2)./max_off_target;
    tv_plot = plot_matrix(:,3)./max_tvs;
    xvec = plot_matrix(:,4);
    figure(1)
    hold on
    plot( xvec, off_t_plot, 'b-', 'LineWidth', 2.0);
    plot( xvec, tv_plot, 'r-', 'LineWidth', 2.0);
    plot( xvec, off_t_plot, 'g.', 'LineWidth', 2.0);
    legend('Off Target Homology','Transcript Variants Matched','Probe Location','location','eastoutside');
    title( [ h_probe_data(1).gene, ' Probe Homology Graph']);
    ylabel('[%]');
    xlabel('Position Along Gene, 5'' to 3'' [BP]');
    axis([0,max(xvec),0,1.1]);
    set(gca, 'FontWeight', 'Bold');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 2.0);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.6, 0.8, 0.4]);
    saveas(gcf,[savepath,'\Homology_',gene_name,'_line'],'epsc');


    % Now go through the matrix and make a list that can be uploaded into blast
    % to double check the results.
    filename = [savepath,'\',gene_name,'_NCBI_Check_File.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tProgram "BLAST210" Output \n\n');
    fprintf(fileID, ['Location\t\t\tmRNA Sequence\n']);
    fprintf(fileID, ['Location\n0 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 0
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n1 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 1
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n2 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 2
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n3 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 3
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fclose(fileID);    
    
    % This makes another list based on the rankings above
    filename = [savepath,'\',gene_name,'_Best_Ranked_Probes.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tProgram "BLAST210" Output: Probes Ranked by Composite Hit Score\n\n');
    fprintf(fileID, ['Location\t\t\tmRNA Sequence\n\n']);

    qst = size(resort_mat,1);
    % Please change this number from 50 to 40 (qst) if you want to obtain results (e.g.,
    % CC16)
    if qst < 50
        numtoprint = resort_mat;
    end
    if qst < numtoprint
        numtoprint = qst;
    end
    for i = 1:numtoprint
        fprintf(fileID, ['>',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, [h_probe_data( resort_mat(i,2) ).mRNA_sequence,'\n']);
    end
    fclose(fileID);
    
    % At this point, nomatter how we got here, it would be good to make a
    % master output sheet which has all the data on it, ranked by the matrix,resort_matrix
    % Reminder: resort_mat = [probesites',probeindicies',probe_hits_sum'];
    filename = [savepath,'\',gene_name,'_Master_Probe_Sheet.txt'];
    fileID = fopen(filename,'w', 'ieee-le','UTF-8');
    % Header
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tMain Output File \n\t\t\t\tAdam M. Andruska, 10/2018\n\n');
    fprintf(fileID, ['\t\t\t\tGene = ',gene_name,'\n']);
    fprintf(fileID, ['\t\t\t\tBitscore Cutoff = ',num2str(bitscore_sig, '%d'),'\n']);
    for i = 1:size(tv_cell,1)
        if i == 1
            fprintf(fileID, ['\t\t\t\tTranscript Variants Considered = \t',tv_cell{i,1},'\n']);
        else
            fprintf(fileID, ['\t\t\t\t\t\t\t\t\t',tv_cell{i,1},'\n']);
        end
    end
    fprintf(fileID,['\n\n\n']);
    % Now cycle and print stuff
    % Please change this number from 50 to 40 (qst) if you want to obtain results (e.g.,
    % CC16)
    qst = size(resort_mat,1);
    if qst < 50
        numtoprint = resort_mat;
    end
    if qst < numtoprint
        numtoprint = qst;
    end
    fprintf(fileID, ['FASTA Sequence IDs and mRNA Sequences For NCBI Blast:\n\n']);
    for i = 1:numtoprint
        fprintf(fileID, ['>',gene_name,'-',num2str( h_probe_data( resort_mat(i,2) ).design_mode,'%d'),'-',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, [h_probe_data( resort_mat(i,2) ).mRNA_sequence,'\n']);
    end

    fprintf(fileID, '\n\n\n');
    for i = 1:numtoprint
        fprintf(fileID, ['\n\tID:  >',gene_name,'-',num2str( h_probe_data( resort_mat(i,2) ).design_mode,'%d'),'-',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, '\tRank\tBP\tmRNA Sequence\t\t\t\t\tBLAST Score\t%%Transcript Variants Targeted\n');
        fprintf(fileID, ['\t',num2str(i,'%d')]);
        fprintf(fileID, ['\t',num2str(resort_mat(i,1),'%d')]);
        fprintf(fileID, ['\t',h_probe_data( resort_mat(i,2) ).mRNA_sequence]);
        fprintf(fileID, ['\t',num2str( resort_mat(i,3),'%d')]);
        fprintf(fileID, ['\t\t',num2str( (h_probe_data( resort_mat(i,2) ).tvs_targeted/size(h_probe_data( resort_mat(i,2) ).total_tvs,1))*100,'%2.1f'),'%%']);

        % This part gives the blast results
        fprintf(fileID, '\n\n\n');
        hits_to_print = 20;
        if hits_to_print > size(h_probe_data(resort_mat(i,2)).ncbi.Hits,2)
            hits_to_print = size(h_probe_data(resort_mat(i,2)).ncbi.Hits,2);
        end
        for j = 1:hits_to_print
            if h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Score >= bitscore_sig
                fprintf(fileID, ['\t\tHit #: ',num2str(j, '%d'),'\n']);
                fprintf(fileID, ['\t\tTarget: ',h_probe_data(resort_mat(i,2)).ncbi.Hits(j).Name,'\n']);
                fprintf(fileID, ['\t\t\tStrand:\t\t',h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Strand,'\n']);
                fprintf(fileID, ['\t\t\tScore:\t\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Score, '%7.4f'),'\n']);
                fprintf(fileID, ['\t\t\tExpect:\t\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Expect, '%6.5e'),'\n']);
                fprintf(fileID, ['\t\t\t# BP Matching:\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Identities.Match, '%d'),'\n']);
                fprintf(fileID, ['\t\t\t%% BP Matching:\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Identities.Percent, '%5.2f'),'\n']);
                fprintf(fileID, ['\t\t\t\tQuery: ',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).QueryIndices(1,1), '%d'),'\t - ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(1,:),' - ',num2str( h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).QueryIndices(1,2), '%d'), '\n']);
                fprintf(fileID, ['\t\t\t\t        \t   ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(2,:), '\n']);
                fprintf(fileID, ['\t\t\t\tSubject: ',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).SubjectIndices(1,1), '%d'),'\t - ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(3,:),' - ', num2str( h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).SubjectIndices(1,2), '%d'), '\n']);
            end
        end
        % This part should print out all the left probes
        fprintf(fileID, ['\n\n\t\t\t\tID\tProbe Name\t\t\tProbe Sequence\t\t\t\t\t\tTm\t\t\t^G\t\t\tHairpins (',num2str(dimer_length,'%d'), ' BP)\t\tDimers (',num2str(dimer_length,'%d'),' BP)']);
        for j = 1:size(h_probe_data( resort_mat(i,2) ).final_h_probes,2)
            if j <= 24
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID,vbseqs(j).name);
                fprintf(fileID, ['\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Left_Name]);
                fprintf(fileID, ['\t\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Full_Left_Probe]);
                fprintf(fileID, ['\t\tTm Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Tm_Left, '%.2f')]);
                fprintf(fileID, ['\t^G Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Delta_G_Left, '%.2f')]);
                fprintf(fileID, ['\tHairpins Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Hairpins_Left, '%d')]);
                fprintf(fileID, ['\tDimers Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Dimers_Left, '%d')]);
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID, ['\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Right_Name]);
                fprintf(fileID, ['\t\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Full_Right_Probe]);
                fprintf(fileID, ['\t\tTm Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Tm_Right, '%.2f')]);
                fprintf(fileID, ['\t^G Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Delta_G_Right, '%.2f')]);
                fprintf(fileID, ['\tHairpins Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Hairpins_Right, '%d')]);
                fprintf(fileID, ['\tDimers Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Dimers_Right, '%d')]);
            else
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID,vbseqs(j).name);
                fprintf(fileID, ['\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Left_Name]);
                fprintf(fileID, ['\t\t\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Full_Left_Probe]);
                fprintf(fileID, ['\t\tTm Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Tm_Left, '%.2f')]);
                fprintf(fileID, ['\t^G Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Delta_G_Left, '%.2f')]);
                fprintf(fileID, ['\tHairpins Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Hairpins_Left, '%d')]);
                fprintf(fileID, ['\tDimers Left  = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Dimers_Left, '%d')]);
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID, ['\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Right_Name]);
                fprintf(fileID, ['\t\t\t',h_probe_data(resort_mat(i,2)).final_h_probes(j).Full_Right_Probe]);
                fprintf(fileID, ['\t\tTm Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Tm_Right, '%.2f')]);
                fprintf(fileID, ['\t^G Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).Delta_G_Right, '%.2f')]);
                fprintf(fileID, ['\tHairpins Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Hairpins_Right, '%d')]);
                fprintf(fileID, ['\tDimers Right = ',num2str(h_probe_data(resort_mat(i,2)).final_h_probes(j).No_Dimers_Right, '%d')]);
            end
        end
    end
    fclose(fileID);
    
elseif load_option == 2
    % This is the option to load pre-existing data:
    gene_name = inputdlg('Enter The Gene Name',...
             'Name Entry', [1 50]);
    gene_name = gene_name{1,1};
    library_directory = 'C:\Users\yalbi\OneDrive\Escritorio\hprobelibrary';
    savepath = [library_directory,'\',gene_name,'_ver_02_2020'];
    S = dir([savepath,'\*.mat']);
    load([savepath,'\',S.name]);
    
    bitscore_sig = 30;
    temp_bit_score = inputdlg('Choose the BitScore ("Score") that defines appropriate homology. A score of 30 or less is the default.',...
                 'Bitscore', [1 50]);
    temp_bit_score = temp_bit_score{1,1};
    if temp_bit_score ~= bitscore_sig
        bitscore_sig = temp_bit_score;
    end
    if ischar(bitscore_sig)
        bitscore_sig = str2num(bitscore_sig);
    end
    
    % Part 0: Ask for an input for all the transcript Variants - we
    % shouldn't have to do this as tvcell is saved.
    
    hit_heatmap = zeros(1, size(h_probe_data,2));
    
    for i = 1:size(h_probe_data,2)
        % At this point the data is reloaded, but the user may have specified a
        % different score to re-analyze the output. Thus, will need to parse
        % all of the probes.
        
        ncbi_short_temp = [];
        for j = 1:size(h_probe_data(i).ncbi.Hits,2)
            
            
            
            name = h_probe_data(i).ncbi.Hits(j).Name;
            strand = h_probe_data(i).ncbi.Hits(j).HSPs.Strand;
            score = h_probe_data(i).ncbi.Hits(j).HSPs.Score;
            expect = h_probe_data(i).ncbi.Hits(j).HSPs.Score;

            tvflag = 0;
            plus_plus_flag = 0;

            % Identify if the hit is important
            if strcmp(strand,'Plus/Plus')

                % Test if the hit is a transcript Variant
                for k = 1:size(tv_cell, 1)
                    if contains(lower(name), tv_cell{k,1})
                        tvflag = 1;
                    end
                end
                plus_plus_flag = 1;
            end

            ncbi_short_temp = [ ncbi_short_temp, [ j; plus_plus_flag; score; tvflag ] ];
        end

        h_probe_data(i).ncbi_short = ncbi_short_temp;   
        h_probe_data(i).ncbi_short_legend = {'Alignment Number'; 'Plus/Plus (1) or Not (0)'; 'Bitscore';'Transcript Variant (1) or not (0)'};


        % Narrow to only Plus/Plus Results
        if isempty(ncbi_short_temp) == 0
            ncbi_report = ncbi_short_temp(:, find(ncbi_short_temp(2,:) == 1));
        else
            off_target_hits = 50;
        end
        if isempty(ncbi_report) == 0
            % Narrow further to only the bitscore
            ncbi_report = ncbi_report(:, find(ncbi_report(3,:) >= bitscore_sig));    
            % Exclude the transcript Variants to get the total number
            % of significant off target hits for the given probe
            off_target_hits = size( ncbi_report(:, find( ncbi_report(4,:) == 0 )), 2);
        else
            off_target_hits = 50;
        end

        % Update the blast reports
        h_probe_data(i).ncbi_report = ncbi_report;
        h_probe_data(i).total_off_target_hits = off_target_hits;

        % Include data for transciprt Variants
        if isempty(ncbi_short_temp) == 0
            h_probe_data(i).tvs_targeted = sum(ncbi_short_temp(4,:));
        else
            h_probe_data(i).tvs_targeted = 0;
        end
        h_probe_data(i).total_tvs = tv_cell;
        h_probe_data(i).significant_bitscore = bitscore_sig;

    end
    
    sorting_matrix = [];
    for i = 1:size(h_probe_data,2)
        sorting_matrix = [sorting_matrix; [i, h_probe_data(i).total_off_target_hits, h_probe_data(i).tvs_targeted, h_probe_data(i).bp_start ] ];
    end
    ranking_matrix = sortrows(sorting_matrix, [2,3],{'ascend','descend'} );
    
    
    % This section will vary the bitscore from a low number to a high
    % number. It will then count up the nubmer of times each probe passes
    % (IE has no significant hits) for a given bitscore. It first has to
    % get rid of any duplicate probesites
    probesites = [-1];
    probeindicies = [-1];
    for i = 1:size(h_probe_data,2)
        tbp = h_probe_data(i).bp_start;
        if ismember(probesites, tbp) == 0
            probesites = [probesites, tbp];
            probeindicies = [probeindicies, i];
        end
    end
    probeindicies = probeindicies(1, 2:size(probeindicies,2));
    probesites = probesites(1, 2:size(probesites,2));
    % Now cycle over bitscore and probesite to find out if the probes are
    % good or not
    probe_hits_sum = zeros(size(probesites));
    for i = bitscore_low:bitscore_high
        for j = 1:size(probesites,2)
            tn = 0;
            temp_mat = h_probe_data(j).ncbi_short;
            if isempty(temp_mat) == 0
                temp_mat = temp_mat(:, find(temp_mat(4,:) == 0) );
                if isempty(temp_mat) == 0
                    temp_mat = temp_mat(3, find(temp_mat(2,:) == 1) );
                    tn = length(find(temp_mat > i));
                end
            end
            if tn > 0
                probe_hits_sum(1,j) =  probe_hits_sum(1,j) + 1;
            end
        end
    end
    resort_mat = [probesites',probeindicies',probe_hits_sum'];
    resort_mat = sortrows(resort_mat,[3],'ascend');
    
    % At this point save the data structure and the mRNA
    filenamemat = [savepath,'\',gene_name,'_DataStructure_',date,'.mat'];
    save(filenamemat, 'h_probe_data','mRNA','ranking_matrix','tv_cell','resort_mat','vbseqs','dimer_length');
    
    
    % USER OUTPUT BELOW
    % This will make a figure showing the probe sites, and respective homology
    plot_matrix = sortrows(ranking_matrix, [4],{'ascend'} );
    max_off_target = max( plot_matrix(:,2) );
    max_tvs = size(h_probe_data(1).total_tvs, 1);
    off_t_plot = plot_matrix(:,2)./max_off_target;
    tv_plot = plot_matrix(:,3)./max_tvs;
    xvec = plot_matrix(:,4);
    figure(1)
    hold on
    plot( xvec, off_t_plot, 'b-', 'LineWidth', 2.0);
    plot( xvec, tv_plot, 'r-', 'LineWidth', 2.0);
    plot( xvec, off_t_plot, 'g.', 'MarkerSize', 8.0);
    legend('Off Target Homology','Transcript Variants Matched','Probe Location','location','eastoutside');
    title( {[ h_probe_data(1).gene, ' Probe Homology Graph'];['BitScore Cutoff = ',num2str(bitscore_sig)]} );
    ylabel('[%]');
    xlabel('Position Along Gene, 5'' to 3'' [BP]');
    axis([0,max(xvec),0,1.1]);
    set(gca, 'FontWeight', 'Bold');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 2.0);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.6, 0.8, 0.4]);
    saveas(gcf,[savepath,'\Homology_',gene_name,'_line'],'epsc');

    figure(2)
    hold on
    bar(probesites, probe_hits_sum, 1.0)
    title( {[ h_probe_data(1).gene, ' Probe Homology Graph'];['BitScore Cutoff = ',num2str(bitscore_sig)]} );
    ylabel('[%]');
    xlabel('Position Along Gene, 5'' to 3'' [BP]');
    axis([0,max(probesites),0,max(probe_hits_sum)]);
    set(gca, 'FontWeight', 'Bold');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 2.0);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.6, 0.8, 0.4]);
    saveas(gcf,[savepath,'\Homology_',gene_name,'_bar'],'epsc');
    

    % Now go through the matrix and make a list that can be uploaded into blast
    % to double check the results.
    filename = [savepath,'\',gene_name,'_NCBI_Check_File.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tProgram "BLAST210" Output \n\n');
    fprintf(fileID, ['Location\t\t\tmRNA Sequence\n']);
    fprintf(fileID, ['Location\n0 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 0
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n1 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 1
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n2 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 2
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fprintf(fileID, ['Location\n3 Off Target Hits\n']);
    for i = 1:size(ranking_matrix, 1)
        if ranking_matrix(i,2) == 3
            fprintf(fileID, ['>',num2str(ranking_matrix(i,4),'%d'),'\n']);
            fprintf(fileID, [h_probe_data( ranking_matrix(i,1) ).mRNA_sequence,'\n']);
        end
    end
    fclose(fileID);
    
    % This makes another list based on the rankings above
    filename = [savepath,'\',gene_name,'_Best_Ranked_Probes.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tProgram "BLAST210" Output: Probes Ranked by Composite Hit Score\n\n');
    fprintf(fileID, ['Location\t\t\tmRNA Sequence\n\n']);

    qst = size(resort_mat,1);
    if qst < 50
        numtoprint = resort_mat;
    end
    if qst < numtoprint
        numtoprint = qst;
    end
    for i = 1:numtoprint
        fprintf(fileID, ['>',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, [h_probe_data( resort_mat(i,2) ).mRNA_sequence,'\n']);
    end
    fclose(fileID);
    
    % At this point, nomatter how we got here, it would be good to make a
    % master output sheet which has all the data on it, ranked by the matrix,resort_matrix
    % Reminder: resort_mat = [probesites',probeindicies',probe_hits_sum'];
    filename = [savepath,'\',gene_name,'_Master_Probe_Sheet.txt'];
    fileID = fopen(filename,'w', 'ieee-le','UTF-8');
    % Header
    fprintf(fileID, '\t\t\t\tPLISH Probe Generator \n');
    fprintf(fileID, '\t\t\t\tMain Output File \n\t\t\t\tAdam M. Andruska, 10/2018\n\n');
    fprintf(fileID, ['\t\t\t\tGene = ',gene_name,'\n']);
    fprintf(fileID, ['\t\t\t\tBitscore Cutoff = ',num2str(bitscore_sig, '%d'),'\n']);
    for i = 1:size(tv_cell,1)
        if i == 1
            fprintf(fileID, ['\t\t\t\tTranscript Variants Considered = \t',tv_cell{i,1},'\n']);
        else
            fprintf(fileID, ['\t\t\t\t\t\t\t\t\t',tv_cell{i,1},'\n']);
        end
    end
    fprintf(fileID,['\n\n\n']);
    % Now cycle and print stuff

    qst = size(resort_mat,1);
    if qst < 50
        numtoprint = resort_mat;
    end
    if qst < numtoprint
        numtoprint = qst;
    end
    fprintf(fileID, ['FASTA Sequence IDs and mRNA Sequences For NCBI Blast:\n\n']);
    for i = 1:numtoprint
        fprintf(fileID, ['>',gene_name,'-',num2str( h_probe_data( resort_mat(i,2) ).design_mode,'%d'),'-',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, [h_probe_data( resort_mat(i,2) ).mRNA_sequence,'\n']);
    end

    fprintf(fileID, '\n\n\n');
    for i = 1:numtoprint
        fprintf(fileID, ['\n\tID:  >',gene_name,'-',num2str( h_probe_data( resort_mat(i,2) ).design_mode,'%d'),'-',num2str(resort_mat(i,1),'%d'),'\n']);
        fprintf(fileID, '\tRank\tBP\tmRNA Sequence\t\t\t\t\tBLAST Score\t%%Transcript Variants Targeted\n');
        fprintf(fileID, ['\t',num2str(i,'%d')]);
        fprintf(fileID, ['\t',num2str(resort_mat(i,1),'%d')]);
        fprintf(fileID, ['\t',h_probe_data( resort_mat(i,2) ).mRNA_sequence]);
        fprintf(fileID, ['\t',num2str( resort_mat(i,3),'%d')]);
        fprintf(fileID, ['\t\t',num2str( (h_probe_data( resort_mat(i,2) ).tvs_targeted/size(h_probe_data( resort_mat(i,2) ).total_tvs,1))*100,'%2.1f'),'%%']);

        % This part gives the blast results
        fprintf(fileID, '\n\n\n');
        hits_to_print = 20;
        if hits_to_print > size(h_probe_data(resort_mat(i,2)).ncbi.Hits,2)
            hits_to_print = size(h_probe_data(resort_mat(i,2)).ncbi.Hits,2);
        end
        for j = 1:hits_to_print
            %cadena = strcmp(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Strand, 'Plus/Plus');
            %if cadena == 0
                if h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Score >= bitscore_sig
            
                    fprintf(fileID, ['\t\tHit #: ',num2str(j, '%d'),'\n']);
                    fprintf(fileID, ['\t\tTarget: ',h_probe_data(resort_mat(i,2)).ncbi.Hits(j).Name,'\n']);
                    fprintf(fileID, ['\t\t\tStrand:\t\t',h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Strand,'\n']);
                    fprintf(fileID, ['\t\t\tScore:\t\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Score, '%7.4f'),'\n']);
                    fprintf(fileID, ['\t\t\tExpect:\t\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Expect, '%6.5e'),'\n']);
                    fprintf(fileID, ['\t\t\t# BP Matching:\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Identities.Match, '%d'),'\n']);
                    fprintf(fileID, ['\t\t\t%% BP Matching:\t',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Identities.Percent, '%5.2f'),'\n']);
                    fprintf(fileID, ['\t\t\t\tQuery: ',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).QueryIndices(1,1), '%d'),'\t - ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(1,:),' - ',num2str( h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).QueryIndices(1,2), '%d'), '\n']);
                    fprintf(fileID, ['\t\t\t\t        \t   ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(2,:), '\n']);
                    fprintf(fileID, ['\t\t\t\tSubject: ',num2str(h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).SubjectIndices(1,1), '%d'),'\t - ', h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).Alignment(3,:),' - ', num2str( h_probe_data(resort_mat(i,2)).ncbi.Hits(j).HSPs(1).SubjectIndices(1,2), '%d'), '\n']);
                end
            %end
        end
        % This part should print out all the left probes
        fprintf(fileID, ['\n\n\t\t\t\tID\tProbe Name\t\t\tProbe Sequence\t\t\t\t\t\tTm\t\t\t^G\t\t\tHairpins (',num2str(dimer_length,'%d'), ' BP)\t\tDimers (',num2str(dimer_length,'%d'),' BP)']);
        for j = 1:size(h_probe_data( resort_mat(i,2) ).final_h_probes,2)
            if j <= 24
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID,vbseqs(j).name);
                fprintf(fileID, ['\t',h_probe_data(i).final_h_probes(j).Left_Name]);
                fprintf(fileID, ['\t\t',h_probe_data(i).final_h_probes(j).Full_Left_Probe]);
                fprintf(fileID, ['\t\tTm Left  = ',num2str(h_probe_data(i).final_h_probes(j).Tm_Left, '%.2f')]);
                fprintf(fileID, ['\t^G Left  = ',num2str(h_probe_data(i).final_h_probes(j).Delta_G_Left, '%.2f')]);
                fprintf(fileID, ['\tHairpins Left  = ',num2str(h_probe_data(i).final_h_probes(j).No_Hairpins_Left, '%d')]);
                fprintf(fileID, ['\tDimers Left  = ',num2str(h_probe_data(i).final_h_probes(j).No_Dimers_Left, '%d')]);
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID, ['\t',h_probe_data(i).final_h_probes(j).Right_Name]);
                fprintf(fileID, ['\t\t',h_probe_data(i).final_h_probes(j).Full_Right_Probe]);
                fprintf(fileID, ['\t\tTm Right = ',num2str(h_probe_data(i).final_h_probes(j).Tm_Right, '%.2f')]);
                fprintf(fileID, ['\t^G Right = ',num2str(h_probe_data(i).final_h_probes(j).Delta_G_Right, '%.2f')]);
                fprintf(fileID, ['\tHairpins Right = ',num2str(h_probe_data(i).final_h_probes(j).No_Hairpins_Right, '%d')]);
                fprintf(fileID, ['\tDimers Right = ',num2str(h_probe_data(i).final_h_probes(j).No_Dimers_Right, '%d')]);
            else
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID,vbseqs(j).name);
                fprintf(fileID, ['\t',h_probe_data(i).final_h_probes(j).Left_Name]);
                fprintf(fileID, ['\t\t\t',h_probe_data(i).final_h_probes(j).Full_Left_Probe]);
                fprintf(fileID, ['\t\tTm Left  = ',num2str(h_probe_data(i).final_h_probes(j).Tm_Left, '%.2f')]);
                fprintf(fileID, ['\t^G Left  = ',num2str(h_probe_data(i).final_h_probes(j).Delta_G_Left, '%.2f')]);
                fprintf(fileID, ['\tHairpins Left  = ',num2str(h_probe_data(i).final_h_probes(j).No_Hairpins_Left, '%d')]);
                fprintf(fileID, ['\tDimers Left  = ',num2str(h_probe_data(i).final_h_probes(j).No_Dimers_Left, '%d')]);
                fprintf(fileID,'\n\t\t\t\t');
                fprintf(fileID, ['\t',h_probe_data(i).final_h_probes(j).Right_Name]);
                fprintf(fileID, ['\t\t\t',h_probe_data(i).final_h_probes(j).Full_Right_Probe]);
                fprintf(fileID, ['\t\tTm Right = ',num2str(h_probe_data(i).final_h_probes(j).Tm_Right, '%.2f')]);
                fprintf(fileID, ['\t^G Right = ',num2str(h_probe_data(i).final_h_probes(j).Delta_G_Right, '%.2f')]);
                fprintf(fileID, ['\tHairpins Right = ',num2str(h_probe_data(i).final_h_probes(j).No_Hairpins_Right, '%d')]);
                fprintf(fileID, ['\tDimers Right = ',num2str(h_probe_data(i).final_h_probes(j).No_Dimers_Right, '%d')]);
            end
        end
    end
    fclose(fileID);
    
    
    % This is the option to make a printout:
else
    % This is the option to load pre-existing data:
    gene_name = inputdlg('Enter The Gene Name',...
             'Name Entry', [1 50]);
    gene_name = gene_name{1,1};
    library_directory = 'C:\Users\yalbi\OneDrive\Documents\GitHub\PLISH_ProbeDesigner\hprobelibrary';
    savepath = [library_directory,'\',gene_name,'_ver_02_2020'];
    S = dir([savepath,'\*.mat']);
    load([savepath,'\',S.name]);

    % User intput for numbers
    candidates = [];
    while isempty(candidates)
        prompts = {'Enter H-Probe Sites (BP) Separated By Spaces:','Enter The Variable Bridges Desired Separated by Spaces. Entering 999 does all bridges.'};
        title = 'Enter the H-Probes You are Intereseted In:';
        dims = [1 50];
        definput = {'0','0'};
        dlg_output = inputdlg(prompts,title,dims,definput);
        candidates = str2num(dlg_output{1,1})';
        bridges_to_print = str2num(dlg_output{2,1})';
    end
    
    if bridges_to_print == 999
        bridges_to_print = (1:28)';
    end

    n_probes = size(h_probe_data,2);
    all_probe_ids = [];
    for i = 1:n_probes
        if h_probe_data(i).design_mode == 2
            all_probe_ids = [all_probe_ids, 999999999];
        else
            all_probe_ids = [all_probe_ids, h_probe_data(i).bp_start];
        end
    end
    final_candidates = [];
    final_probe_i = [];
    for i = 1:size(candidates,1);
        if ismember(candidates(i,1), all_probe_ids)
            final_candidates = [final_candidates, candidates(i,1)];
            final_probe_i = [final_probe_i, find(all_probe_ids == candidates(i,1))];
        end
    end


    % By this point, we have made a proofread list of the probes that the user
    % is intersted in. Now to conjugate them and make a report.

    table_1_top = {     'Count','Fluorophore','Common Circle','Probe Name','H-Probe Sequence','mRNA Target','Top Blast Result (Check)' };

    if strcmp('h',lower(gene_name(1,1)))
        species_pre = 'hum';
    elseif strcmp('m',lower(gene_name(1,1)))
        species_pre = 'mm';
    else
        species_pre = 'xx';
    end

    new_circle_name = 'CCC2.1';
    old_circle_name = 'CC3A.1';
    
    
    for i = 1:size(bridges_to_print,1)
        
        for j = 1:size(final_probe_i,2)

            temp_vb = bridges_to_print(i,1);
            temp_site = final_probe_i(1,j);
            
            if temp_vb < 25
                circle_temp = new_circle_name;
            else
                circle_temp = old_circle_name;
            end

            table_add = {...
                    num2str(j),...
                    vbseqs(temp_vb).name,...
                    circle_temp,...
                    h_probe_data(temp_site).final_h_probes(temp_vb).Left_Name,...
                    h_probe_data(temp_site).final_h_probes(temp_vb).Full_Left_Probe,...
                    h_probe_data(temp_site).mRNA_sequence,...
                    h_probe_data(temp_site).ncbi.Hits(1).Name;...
                    ' ',...
                    ' ',...
                    ' ',...
                    h_probe_data(temp_site).final_h_probes(temp_vb).Right_Name,...
                    h_probe_data(temp_site).final_h_probes(temp_vb).Full_Right_Probe,...
                    ' ',...
                    ' ' };

              table_1_top = [ table_1_top; table_add];

        end
        table_1_top = [ table_1_top; { ' ',' ',' ',' ',' ',' ',' ' } ];
    end

    filename = [savepath,'\',gene_name,'_HProbeSheet_',num2str(final_candidates),'.xlsx'];
    xlswrite(filename,table_1_top);
    h = msgbox(['The output file should be found in ',savepath,'.'], 'Program "Printout" Is Done','help');
    
end

% fin

