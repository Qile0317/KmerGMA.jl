# file for the UPCOMING feature to check for RSS's in the buffer of hits, specifically for immunoglobulins

using BioSequences, BioAlignments

export HumanRSSV
export HumanRSSD
export Align_RSS
export RSS_dist
export is_RSS

rand_12nt_spacer = dna"N"^12
rand_23nt_spacer = dna"N"^23

const HumanRSSV = dna"CACAGTG" * rand_12nt_spacer * dna"ACAAAAACC"
const HumanRSSD = dna"CACAGTG" * rand_23nt_spacer * dna"ACAAAAACC"

# seq should be the last buffer 
function Align_RSS(seq::Seq, RSS::Seq = HumanRSSV; score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1))
    return pairalign(SemiGlobalAlignment(), RSS, seq, score_model)
end

function RSS_dist(RSS1::Seq, RSS2::Seq = HumanRSSV)
    dist = 0; for i in eachindex(RSS1)
        if RSS1[i] != RSS2[i]; dist += 1 end
    end; return dist 
end

# evaluate if sequence is a true RSS. currently just based on hamming but other ways are possible too
function is_RSS(RSS_align::AlignObj, thr::Int = 1)
    RSS_candidate = RSS_align.aln.b[cigar_to_UnitRange(RSS_align)]
    return RSS_dist(RSS_candidate, RSS_align.aln.a.seq) <= thr
end

# funnily enough, KmerGMA can technically also do this