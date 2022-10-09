from collections import Counter


# might split this schema into tmp/final TODO
class Schema:  # TODO ???drop 'pos_in_read', 'read_names'????
    CHROMOSOME = "seq"
    POSITION = "pos"
    REFERENCE = "ref"
    READS = "reads"
    RESULTS = "res"
    QUALITY = "qual"
    POSITION_IN_READ = "pos_in_read"
    READ_NAMES = "read_names"
    FLAGS = "flags"
    MAPPING_QUALITY = "map"
    STRUCTURAL_VARIANTS = "SV"
    RESULTS_NO_CARET = "res_no_caret"
    RESULTS_NUCLEOTIDES = "res_nuc"
    RESULTS_NUCLEOTIDES_FILTERED = "res_nuc_filtered"
    QUALITY_FILTERED = "qual_filtred"
    RESULTS_INDELS = "res_indel"
    SINGLE_NUCLEOTIDE_CALLS = "SNV_calls"
    A_CALLS = "A_calls"
    T_CALLS = "T_calls"
    C_CALLS = "C_calls"
    G_CALLS = "G_calls"
    A_PVALUE = "A_pvalue"
    T_PVALUE = "T_pvalue"
    C_PVALUE = "C_pvalue"
    G_PVALUE = "G_pvalue"
    INDEL_CALLS = "INDEL_calls"
    STRUCTURAL_CALLS = "SV_calls"
    INDEL_PVALUE = "INDEL_pvalue"
    STRUCTURAL_PVALUE = "SV_pvalue"
    HOMOPOLYMER = "In_or_ajacent_to_homopolymer_of_length"
    READ_DEPTH_POST_FILTER = "read_depth_post_filtering"


def convert_to_upper(res: str, ref: str) -> str:
    """converts results str to upper case, replaces "," and "." to reference str
    and removes "*" to allow phred scores to align"""
    assert ref.isupper()
    return res.replace(",", ref).replace(".", ref).replace("*", "").upper()


def remove_ref_from_tumour_res(
    normal_results: str, tumour_results: str, ref: str
) -> str:
    """Removes any reference seq from the tumour results str"""
    most_common_base_in_normal = Counter(
        convert_to_upper(normal_results, ref)
    ).most_common()[0][0]
    return convert_to_upper(tumour_results, ref).replace(most_common_base_in_normal, "")
