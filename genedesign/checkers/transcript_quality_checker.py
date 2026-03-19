from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.seq_utils.hairpin_counter import hairpin_counter


class TranscriptQualityChecker:
    """
    Composite checker used only inside TranscriptDesigner.

    It still uses the existing benchmark checkers, but returns richer metadata
    so TranscriptDesigner can search more intelligently, with extra emphasis
    on identifying and repairing hairpin-prone windows.
    """

    def __init__(self):
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None

        self.chunk_size = 50
        self.step = 25
        self.min_stem = 3
        self.min_loop = 4
        self.max_loop = 9

    def initiate(self):
        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()

        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

    def _hairpin_window_details(self, transcript_dna: str) -> list[dict]:
        """
        Mirror the benchmark hairpin logic, but keep detailed counts so the
        designer can target the worst windows.
        """
        details = []

        if len(transcript_dna) < self.chunk_size:
            return details

        for start in range(0, len(transcript_dna) - self.chunk_size + 1, self.step):
            end = start + self.chunk_size
            chunk = transcript_dna[start:end]
            hairpin_count, hairpin_string = hairpin_counter(
                chunk,
                self.min_stem,
                self.min_loop,
                self.max_loop,
            )
            if hairpin_count > 1:
                details.append(
                    {
                        "start": start,
                        "end": end,
                        "count": hairpin_count,
                        "string": hairpin_string,
                    }
                )

        return details

    def run(self, utr: str, codons: list[str]) -> dict:
        cds = "".join(codons).upper()
        transcript_dna = (utr + cds).upper()

        codons_ok, diversity, rare_count, cai = self.codon_checker.run(codons)
        forbidden_ok, forbidden_hit = self.forbidden_checker.run(transcript_dna)
        promoter_ok, promoter_hit = self.promoter_checker.run(transcript_dna)
        hairpin_ok, hairpin_hit = hairpin_checker(transcript_dna)

        if hairpin_ok:
            hairpin_details = []
        else:
            hairpin_details = self._hairpin_window_details(transcript_dna)

        hairpin_bad_windows = [(d["start"], d["end"]) for d in hairpin_details]
        hairpin_total_excess = sum(max(0, d["count"] - 1) for d in hairpin_details)
        hairpin_max_window_count = max([d["count"] for d in hairpin_details], default=0)

        hard_failures = (
            int(not forbidden_ok)
            + int(not promoter_ok)
            + int(not hairpin_ok)
            + int(not codons_ok)
        )

        # Heavily prioritize fewer hard failures and fewer / weaker hairpin windows.
        rank = (
            hard_failures,
            hairpin_total_excess,
            len(hairpin_bad_windows),
            hairpin_max_window_count,
            0 if forbidden_ok else 1,
            0 if promoter_ok else 1,
            0 if codons_ok else 1,
            rare_count,
            -cai,
            -diversity,
        )

        score = 0.0
        score -= 1500.0 * hard_failures
        score -= 400.0 * hairpin_total_excess
        score -= 150.0 * len(hairpin_bad_windows)
        score -= 50.0 * max(0, hairpin_max_window_count - 1)

        score += 125.0 if forbidden_ok else -250.0
        score += 125.0 if promoter_ok else -250.0
        score += 125.0 if hairpin_ok else -250.0
        score += 80.0 if codons_ok else -120.0

        score += 220.0 * cai
        score += 120.0 * diversity
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
            "hairpin_bad_window_details": hairpin_details,
            "hairpin_total_excess": hairpin_total_excess,
            "hairpin_max_window_count": hairpin_max_window_count,
        }
