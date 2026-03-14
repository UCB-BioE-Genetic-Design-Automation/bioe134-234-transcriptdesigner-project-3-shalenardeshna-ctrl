from genedesign.seq_utils.hairpin_counter import hairpin_counter


def hairpin_checker(dna):
    """
    Checks for bad hairpin structures in the DNA sequence.

    Improvements over the original version:
    - sequences shorter than 50 bp are still checked
    - the tail/end of long sequences is always checked
    - returns False if any checked window has >1 hairpin
    """
    dna = dna.upper()

    chunk_size = 50
    step = 25
    min_stem = 3
    min_loop = 4
    max_loop = 9

    if len(dna) == 0:
        return True, None

    # Short sequences should still be checked.
    if len(dna) <= chunk_size:
        hairpin_count, hairpin_string = hairpin_counter(dna, min_stem, min_loop, max_loop)
        if hairpin_count > 1:
            return False, hairpin_string
        return True, None

    # Check overlapping windows across the sequence.
    starts = list(range(0, len(dna) - chunk_size + 1, step))

    # Always include the final window so the tail is not skipped.
    final_start = len(dna) - chunk_size
    if starts[-1] != final_start:
        starts.append(final_start)

    for i in starts:
        chunk = dna[i:i + chunk_size]
        hairpin_count, hairpin_string = hairpin_counter(chunk, min_stem, min_loop, max_loop)

        if hairpin_count > 1:
            return False, hairpin_string

    return True, None


if __name__ == "__main__":
    result, hairpin = hairpin_checker(
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCAAAAAAAGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
    print(result, hairpin)
