using BioAlignments

const AlignResult = PairwiseAlignmentResult{Int64, LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}
const SubAlignResult = PairwiseAlignmentResult{Int64,LongSequence{DNAAlphabet{4}},LongSubSeq{DNAAlphabet{4}}}
const AlignObj = Union{AlignResult, SubAlignResult}

const char_ints = Set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])

@inline function cigar_to_UnitRange(aligned_obj)
    cigar_str = cigar(aligned_obj.aln.a.aln)
    curr_num, char_count, num_sum, lower = 0, 0, 0, 0
    n = length(cigar_str)
    for i in 1:n
        if i == n; return (lower + 1):num_sum; end 
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

# helper to get the aligned hit unitrange depending 
@inline function align_unitrange(
    seq::DnaSeq, seq_UnitRange::UnitRange{Int},
    consensus_seq::DnaSeq, windowsize::Int,
    sequence_length::Int,
    score_model::AffineGapScoreModel{Int} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1))

    aligned_obj = pairalign(SemiGlobalAlignment(),
                        view(consensus_seq, 1:windowsize),
                        view(seq, seq_UnitRange),
                        score_model)

    aligned_UnitRange = cigar_to_UnitRange(aligned_obj)
    return max(
        1, first(seq_UnitRange) + first(aligned_UnitRange) - 1):min(
            first(seq_UnitRange) + last(aligned_UnitRange) - 1, sequence_length)
end

export align_unitrange

# known issue in alignment: sometimes the C/G in a V gene doesn match to the correct corresponding G/C, messing up the alignment and cigar