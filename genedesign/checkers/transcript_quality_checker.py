from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker


class TranscriptQualityChecker:
    """
    Composite checker used only inside TranscriptDesigner.

    It still uses the existing benchmark checkers, but returns richer metadata
    so TranscriptDesigner can search more intelligently.
    """

    def __init__(self):
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None

    def initiate(self):
        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()

        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

    def _hairpin_window_failures(self, transcript_dna: str, chunk_size: int = 50, step: int = 25):
        """
        Count failing windows using the existing hairpin_checker.
        """
        bad_windows = []

        if len(transcript_dna) <= chunk_size:
            passed, _ = hairpin_checker(transcript_dna)
            if not passed:
                bad_windows.append((0, len(transcript_dna)))
            return bad_windows

        for start in range(0, len(transcript_dna) - chunk_size + 1, step):
            end = start + chunk_size
            chunk = transcript_dna[start:end]
            passed, _ = hairpin_checker(chunk)
            if not passed:
                bad_windows.append((start, end))

        return bad_windows

    def run(self, utr: str, codons: list[str]) -> dict:
        """
        Parameters
        ----------
        utr : str
            RBS UTR sequence.
        codons : list[str]
            CDS codons including stop codon.

        Returns
        -------
        dict
            Summary of all checks plus ranking data for search.
        """
        cds = "".join(codons).upper()
        transcript_dna = (utr + cds).upper()

        codons_ok, diversity, rare_count, cai = self.codon_checker.run(codons)
        forbidden_ok, forbidden_hit = self.forbidden_checker.run(transcript_dna)
        promoter_ok, promoter_hit = self.promoter_checker.run(transcript_dna)
        hairpin_ok, hairpin_hit = hairpin_checker(transcript_dna)

        hairpin_bad_windows = self._hairpin_window_failures(transcript_dna)
        hard_failures = (
            int(not forbidden_ok)
            + int(not promoter_ok)
            + int(not hairpin_ok)
            + int(not codons_ok)
        )

        # Lower hard-failure count dominates.
        # Then reduce number of bad hairpin windows.
        # Then improve codon metrics.
        rank = (
            hard_failures,
            0 if hairpin_ok else 1,
            len(hairpin_bad_windows),
            0 if forbidden_ok else 1,
            0 if promoter_ok else 1,
            0 if codons_ok else 1,
            rare_count,
            -diversity,
            -cai,
        )

        score = 0.0
        score -= 1000.0 * hard_failures
        score -= 150.0 * len(hairpin_bad_windows)

        score += 100.0 if forbidden_ok else -200.0
        score += 100.0 if promoter_ok else -200.0
        score += 100.0 if hairpin_ok else -200.0
        score += 100.0 if codons_ok else -200.0

        score += 250.0 * diversity
        score += 300.0 * cai
        score -= 40.0 * rare_count

        passed = forbidden_ok and promoter_ok and hairpin_ok and codons_ok

        return {
            "passed": passed,
            "score": score,
            "rank": rank,
            "hard_failures": hard_failures,
            "codons_ok": codons_ok,
            "diversity": diversity,
            "rare_count": rare_count,
            "cai": cai,
            "forbidden_ok": forbidden_ok,
            "forbidden_hit": forbidden_hit,
            "promoter_ok": promoter_ok,
            "promoter_hit": promoter_hit,
            "hairpin_ok": hairpin_ok,
            "hairpin_hit": hairpin_hit,
            "hairpin_bad_windows": hairpin_bad_windows,
        }
