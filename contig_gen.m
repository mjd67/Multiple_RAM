function [contigs] = contig_gen(Fwd,Rev,set_length)
%Contig generation function, 9/3/25
% %**throws out lengths less than 235 for now
%There seems to be some read error correction I need to do...

contigs = strings(set_length,1);
for j =1:set_length
    top=Fwd(j).Sequence;
    bottom = Rev(j).Sequence;
    
    if(length(top) >235 && length(bottom) > 235)
        %Take in the forward (top) and reverse (bottom) reads and generates a
        %single contig. Bottom comes in as reverse compliment  
        %Overlap of k-mer length 8 at the 225-233 lengths
        bottom = seqrcomplement(bottom);
        k=strfind(bottom,top(225:233));
        contigs(j)=strcat(top(1:224),bottom(k:length(bottom)));
    end

end

contigs(contigs == "") = [];


end