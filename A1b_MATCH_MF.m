function match = MATCH_MF(ml,snp)
%Looks for entries in masterlist "ml" that are also present in the IndvMutList of the replicate "snp".
%
%ml and snp must be structures with the fields pos (the position of the
%SNP) and atcg (the alternate bp for that position)
%
%OUTPUT: match contains the positions in the genome that are identical
%(position and alternate bp) in both structures.

[match_pos, match_posidx]=ismember(ml.pos, snp.pos);        % With positions are on both lists ? match_pos is a mask (logicals) for the masterlist, match_posidx gives the index on the snplist of the replicate
match_nt = strcmp(ml.atcg(match_pos), snp.atcg(nonzeros(match_posidx)));    % Are the detected alternate alleles the same as the expected ? match_nt gives a mask for the 

paren = @(x, varargin) x(varargin{:});
match = paren(ml.pos(match_pos), match_nt); 
