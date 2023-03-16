using BioAlignments

const char_ints = Set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])

function cigar_to_UnitRange(aligned_obj)
    cigar_str = cigar(aligned_obj.aln.a.aln)
    curr_num, char_count, num_sum, lower = 0, 0, 0, 0
    n = length(cigar_str)
    for i in 1:n
        if i == n; return (lower+1):num_sum; end 
        if cigar_str[i] in char_ints
            curr_num = (curr_num * 10) + parse(Int, cigar_str[i])
        else
            char_count += 1
            if char_count == 1; lower = curr_num end
            num_sum += curr_num
            curr_num = 0
        end
    end
end

export cigar_to_UnitRange

# known issue in alignment: sometimes the C/G in a V gene doesn match to the correct corresponding G/C, messing up the alignment and cigar