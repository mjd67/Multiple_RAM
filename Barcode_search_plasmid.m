function exp = Barcode_search_plasmid(FASTQStructR,bar_seq_DNA,full_run_count)
%inputs

%** needs to be updated for non-e coli, removed that section for now

% FASTQStructR
% bar_seq_DNA
% full_run_count
% is the seq run from RNA (true or false)


%FASTQStructR are the structure of the reverse reads
%excel_data is the list of barcodes
%full_run_count is the total number of reads in the seq run
%s is how large of the catch all count should be, default to 10000 to be
%safe

pat = 'AGTGCAAAGATGACGGG'; %Stem, seen
DNA_pat  = 'TGCGAAAGTATATTGATTAG'; %Ribozyme sequence
%exp 1-10 =  order from the excel file
% good = 0; %exp 11
% empty = 0; %exp 12
% short = 0; %exp 13
% chimeric=0; %exp 14
% misscall = 0;%exp 15
FASTQStructR = fastqread(FASTQStructR);
exp = zeros(1,15);
for i = 1:full_run_count
    seq = seqrcomplement(FASTQStructR(1,i).Sequence); %looking at the reverse reads for barcodes, so swapping orientation
    k = strfind(seq,pat);
    if isempty(k)%exact match, maybe can give 1 mismatch later
        exp(12) = exp(12) +1;
    elseif k+27 > length(seq) %this may be wrong
        exp(13) = exp(13) +1;
    elseif ~contains(seq, DNA_pat)  %checking for a chimeric sequence, RNA only
        exp(14) =  exp(14) +1;
    else
        barcode = seq(k+17:k+26);
        look = find(strcmp(barcode, bar_seq_DNA));
        if isempty(look)
            exp(15) = exp(15)+1;
        else
            exp(11) = exp(11)+1;
            exp(look) = exp(look) + 1;
        end
    end
    
end
