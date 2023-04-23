using BioAlignments

const AlignResult = PairwiseAlignmentResult{Int64, LongSubSeq{DNAAlphabet{4}}, LongSubSeq{DNAAlphabet{4}}}
const RSSAlignResult = PairwiseAlignmentResult{Int64,LongSequence{DNAAlphabet{4}}, LongSubSeq{DNAAlphabet{4}}}

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

# helper to get the aligned hit unitrange. returns the new seq_UnitRange
@inline function align_unitrange(
    seq::DnaSeq, seq_UnitRange::UnitRange{Int},
    consensus_seq::DnaSeq, windowsize::Int,
    sequence_length::Int,
    score_model::AffineGapScoreModel{Int} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    do_return_align::Bool = false,
    result_align_vec::Vector{AlignResult} = AlignResult[])

    aligned_obj = pairalign(SemiGlobalAlignment(),
                        view(consensus_seq, 1:windowsize),
                        view(seq, seq_UnitRange),
                        score_model)

    if do_return_align; push!(result_align_vec, aligned_obj) end

    aligned_UnitRange = cigar_to_UnitRange(aligned_obj)
    return max(
        1, first(seq_UnitRange) + first(aligned_UnitRange) - 1):min(
            first(seq_UnitRange) + last(aligned_UnitRange) - 1, sequence_length)
end

export align_unitrange

# known issue in alignment: sometimes the C/G in a V gene doesn match to the correct corresponding G/C, messing up the alignment and cigar

# convinience function to append a hit - slows down the main function slightly but its speed is negligible
@inline function append_hit!(
    resultVec::Vector{FASTA.Record},
    record::FASTA.Record,
    seq::Seq,
    do_overlap::Bool,
    windowsize::Int,
    currminim::Float64,
    seq_UnitRange::UnitRange{Int},
    genome_pos::Int)

    # uses teeny bit more memory to reduce code repetition
    MatchPos_UnitRange::UnitRange{Int} = seq_UnitRange
    if do_overlap
        MatchPos_UnitRange = (first(seq_UnitRange) - windowsize + 1):(last(seq_UnitRange) - windowsize + 1)
    end
    push!(
        resultVec, 
        FASTA.Record(
            FASTA.identifier(record) *
                " | dist = " * string(round(currminim, digits = 2)) *
                " | MatchPos = $MatchPos_UnitRange" *
                " | GenomePos = $genome_pos"*
                " | Len = "*string(last(seq_UnitRange)-first(seq_UnitRange) + 1), 
            view(seq, seq_UnitRange)
    ))
end

export append_hit!