function exp = Barcode_search_community(FASTQStructR,FASTQStructF,bar_seq_DNA,full_run_count)
    %inputs
    
    %** needs to be updated for non-e coli, removed that section for now
    
    % FASTQStructR
    % bar_seq_DNA
    % full_run_count
    % is the seq run from RNA (true or false)
    
    
    %FASTQStructR are the structure of the reverse reads, file is generated
    %here in the script. Will also need to forward read for species
    %identification
    %excel_data is the list of barcodes
    %full_run_count is the total number of reads in the seq run
    %s is how large of the catch all count should be, default to 10000 to be
    %safe
    
    pat = 'AGTGCAAAGATGACGGG'; %Stem, seen
    DNA_pat  = 'TGCGAAAGTATATTGATTAG'; %Ribozyme sequence
    
    
    e_coli16S  = "GCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCA";
    vmax_16S_1  = upper("gcaacgcgaagaaccttacctactcttgacatccagagaacttttcagagatgaattggtgccttcgggaactctgagacaggtgctgcatggctgtcgtcagctcgtgttgtgaaatgttgggttaagtcccgcaacgagcgcaacccttatccttgtttgccagcgagtaatgtcgggaactccagggagactgccggtgataaaccggaggaaggtggggacgacgtcaagtcatcatggcccttacgagtagggctacacacgtgctacaatggcgcatacagagggcggccaacttgcgaaagtgagcgaatcccaaaaagtgcgtcgtagtccggattggagtctgcaactcgactccatgaagtcggaatcgctagtaatcgtggatcagaatgccacggtgaatacgttcccgggccttgtacacaccgcccgtca");
    %vmax_16S_2  = upper("gcaacgcgaagaaccttacctactcttgacatccagagaactttccagagatggattggtgccttcgggaactctgagacaggtgctgcatggctgtcgtcagctcgtgttgtgaaatgttgggttaagtcccgcaacgagcgcaacccttatccttgtttgccagcgagtaatgtcgggaactccagggagactgccggtgataaaccggaggaaggtggggacgacgtcaagtcatcatggcccttacgagtagggctacacacgtgctacaatggcgcatacagagggcggccaacttgcgaaagtgagcgaatcccaaaaagtgcgtcgtagtccggattggagtctgcaactcgactccatgaagtcggaatcgctagtaatcgtggatcagaatgccacggtgaatacgttcccgggccttgtacacaccgcccgtca");
    %vmax_16S_3  = upper("gcaacgcgaagaaccttacctactcttgacatccagagaattttccagagatggattagtgccttcgggaactctgagacaggtgctgcatggctgtcgtcagctcgtgttgtgaaatgttgggttaagtcccgcaacgagcgcaacccttatccttgtttgccagcgagtaatgtcgggaactccagggagactgccggtgataaaccggaggaaggtggggacgacgtcaagtcatcatggcccttacgagtagggctacacacgtgctacaatggcgcatacagagggcggccaacttgcgaaagtgagcgaatcccaaaaagtgcgtcgtagtccggattggagtctgcaactcgactccatgaagtcggaatcgctagtaatcgtggatcagaatgccacggtgaatacgttcccgggccttgtacacaccgcccgtca");
    putida_16S = "GCAACGCGAAGAACCTTACCAGGCCTTGACATGCAGAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTCTGACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTACCAGCACGTTATGGTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGGCCTGGGCTACACACGTGCTACAATGGTCGGTACAGAGGGTTGCCAAGCCGCGAGGTGGAGCTAATCTCACAAAACCGATCGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGAATCAGAATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACC";
    
    S16_seqs = [e_coli16S,vmax_16S_1,putida_16S];
    
    %exp 1-10 =  order from the excel file
    % good = 0; %exp 11
    % empty = 0; %exp 12
    % short = 0; %exp 13
    % chimeric=0; %exp 14
    % misscall = 0;%exp 15
    % coli %exp16
    % vibrio %exp17
    % putida %exp18
    % removed for species identity exp 19;
    % coli 21-30 ; vibrio 31-40 ; putida 41-50
    FASTQStructR = fastqread(FASTQStructR);
    FASTQStructF = fastqread(FASTQStructF);
    exp = zeros(1,50);
    for i = 1:full_run_count
        species = false;
        seq = seqrcomplement(FASTQStructR(1,i).Sequence); %looking at the reverse reads for barcodes, so swapping orientation
        k = strfind(seq,pat);
        if isempty(k)%exact match, maybe can give 1 mismatch later
            exp(12) = exp(12) +1;
        elseif k+27 > length(seq) %this may be wrong
            exp(13) = exp(13) +1;
        elseif contains(seq, DNA_pat)  %checking for a chimeric sequence, RNA only
            exp(14) =  exp(14) +1;
        else
            barcode = seq(k+17:k+26);
            look = find(strcmp(barcode, bar_seq_DNA));
            if isempty(look)
                exp(15) = exp(15)+1;
            else
                exp(11) = exp(11)+1;
                exp(look) = exp(look) + 1;
                species = true;
            end
        end
    
    
        %species matching
        if species == true
            seqR = FASTQStructR(1,i).Sequence;
            seqF = FASTQStructF(1,i).Sequence;
            a = length(seqF) < 225 ; %removing reads that are too short
            b = length(seqR) < 225; %removing reads that are too short
            c= isempty(strfind(seqF,upper('gcaacgcgaagaaccttacc')));%fwd primer
            d = isempty(strfind(seqR,'ATTTGAACCGACGATCTTCGG')); %rev primer
            e = ~contains(seqrcomplement(seqR),'AGTGCAAAGATGACGGG'); %stem
            if  a || b || c || d || e
                exp(19) = exp(19) +1; %maybe look into breaking down each of these if cell 19 is more than 10% of reads?
            else
                k = strfind(seqF,upper('gcaacgcgaagaaccttacc')) +20;
                r = strfind(seqR,'ATTTGAACCGACGATCTTCGGG')+106; %switching from 107-106 to maybe help vibrio?
                
                top = seqF(k:end);
                bottom = seqrcomplement(seqR(r:end)); 
        
                %space  = 12 + (230 - length(top)) + (143-length(bottom));
                nnn=strcat(repmat('N',1,12));
                ftigs = {strcat(top,nnn,bottom)};
    
                scores = zeros(3,1);
                for m = 1:3
                    match= char(S16_seqs(m));
                    scores(m) = swalign(ftigs{1},match(20:(end-2)));
                end
                [~, idx] = max(scores);
                dilate = (idx+1)*10 + look;
                exp(dilate) = exp(dilate) +1;
    
            end
        end
    end

end