from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker


class TranscriptQualityChecker:
    """
    Composite checker for a transcript design.

    Evaluates:
    - codon quality (diversity, rare codons, CAI)
    - forbidden sequences
    - internal promoters
    - hairpins

    This checker is meant to be used inside TranscriptDesigner to compare
    candidate designs. It does NOT replace the benchmark's own checks.
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
            Summary of all checks plus an overall score.
        """
        cds = "".join(codons)
        transcript_dna = (utr + cds).upper()

        # CodonChecker expects a list of codons.
        codons_ok, diversity, rare_count, cai = self.codon_checker.run(codons)

        forbidden_ok, forbidden_hit = self.forbidden_checker.run(transcript_dna)
        promoter_ok, promoter_hit = self.promoter_checker.run(transcript_dna)
        hairpin_ok, hairpin_hit = hairpin_checker(transcript_dna)

        # Large penalties for hard failures, rewards for good coding properties.
        score = 0.0
        score += 100.0 if forbidden_ok else -200.0
        score += 100.0 if promoter_ok else -200.0
        score += 60.0 if hairpin_ok else -120.0
        score += 25.0 if codons_ok else -25.0
        score += 40.0 * cai
        score += 10.0 * diversity
        score -= 8.0 * rare_count

        passed = forbidden_ok and promoter_ok and hairpin_ok and codons_ok

        return {
            "passed": passed,
            "score": score,
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
        }
