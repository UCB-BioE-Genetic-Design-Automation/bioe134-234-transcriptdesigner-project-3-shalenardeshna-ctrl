from collections import defaultdict

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.transcript_quality_checker import TranscriptQualityChecker
from genedesign.checkers.hairpin_checker import hairpin_checker


class TranscriptDesigner:
    """
    Reverse-translates a protein sequence into a DNA sequence and chooses
    an RBS while trying to satisfy quality constraints.

    This implementation uses a sliding-window local search over synonymous
    codons instead of a single small global greedy pass.
    """

    def __init__(self):
        self.codon_choices = {}
        self.rbsChooser = None
        self.qualityChecker = None

        # Sliding-window search settings
        self.window_codons = 12
        self.window_step = 4
        self.max_passes = 4
        self.local_flank_nt = 72
        self.max_branch = 4
        self.stop_codon = "TAA"

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.qualityChecker = TranscriptQualityChecker()
        self.qualityChecker.initiate()

        self._init_codon_choices()

    def _init_codon_choices(self) -> None:
        """
        Build synonymous codon lists ordered by codon usage preference while
        filtering rare codons when possible.
        """
        genetic_code = {
            "A": ["GCT", "GCC", "GCA", "GCG"],
            "C": ["TGT", "TGC"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "F": ["TTT", "TTC"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "H": ["CAT", "CAC"],
            "I": ["ATT", "ATC", "ATA"],
            "K": ["AAA", "AAG"],
            "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            "M": ["ATG"],
            "N": ["AAT", "AAC"],
            "P": ["CCT", "CCC", "CCA", "CCG"],
            "Q": ["CAA", "CAG"],
            "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            "T": ["ACT", "ACC", "ACA", "ACG"],
            "V": ["GTT", "GTC", "GTA", "GTG"],
            "W": ["TGG"],
            "Y": ["TAT", "TAC"],
        }

        codon_checker = self.qualityChecker.codon_checker
        freqs = codon_checker.codon_frequencies
        rare_threshold = codon_checker.rare_codon_threshold

        self.codon_choices = {}
        for aa, codons in genetic_code.items():
            usable = [codon for codon in codons if freqs.get(codon, 0.0) >= rare_threshold]
            if not usable:
                usable = codons[:]
            usable.sort(key=lambda codon: freqs.get(codon, 0.0), reverse=True)
            self.codon_choices[aa] = usable

    def _available_rbs_options(self, ignores: set):
        options = [opt for opt in self.rbsChooser.rbsOptions if opt not in ignores]
        if not options:
            raise Exception("No valid RBS option available.")
        return options

    def _seed_codons(self, peptide: str) -> list[str]:
        """
        Create a high-CAI but still somewhat diverse initial CDS by rotating
        through the best few synonymous codons for repeated amino acids.
        """
        counts = defaultdict(int)
        codons = []

        for aa in peptide:
            choices = self.codon_choices[aa]
            rotate = min(3, len(choices))
            choice_index = counts[aa] % rotate
            codons.append(choices[choice_index])
            counts[aa] += 1

        return codons

    def _window_ranges(self, n_positions: int) -> list[tuple[int, int]]:
        """
        Return overlapping sliding windows over codon positions.
        """
        if n_positions <= self.window_codons:
            return [(0, n_positions)]

        windows = []
        start = 0

        while True:
            end = min(n_positions, start + self.window_codons)
            if end - start < self.window_codons:
                start = max(0, n_positions - self.window_codons)
                end = n_positions

            if windows and windows[-1] == (start, end):
                break

            windows.append((start, end))

            if end >= n_positions:
                break

            start += self.window_step

        return windows

    def _transcript_dna(self, utr: str, codons: list[str]) -> str:
        return (utr + "".join(codons)).upper()

    def _codon_score(self, codons: list[str]) -> float:
        codons_ok, diversity, rare_count, cai = self.qualityChecker.codon_checker.run(codons)
        score = 0.0
        score += 25.0 if codons_ok else -25.0
        score += 40.0 * cai
        score += 10.0 * diversity
        score -= 8.0 * rare_count
        return score

    def _local_structural_score(self, utr: str, codons: list[str], pos: int) -> tuple[float, bool]:
        """
        Score only the local region around one codon change. This keeps the
        sliding-window search much faster than rescoring the full transcript
        after every single synonymous substitution.
        """
        dna = self._transcript_dna(utr, codons)
        center = len(utr) + (pos * 3) + 1
        start = max(0, center - self.local_flank_nt)
        end = min(len(dna), center + self.local_flank_nt)
        segment = dna[start:end]

        forbidden_ok, _ = self.qualityChecker.forbidden_checker.run(segment)
        promoter_ok, _ = self.qualityChecker.promoter_checker.run(segment)
        hairpin_ok, _ = hairpin_checker(segment)

        score = 0.0
        score += 100.0 if forbidden_ok else -200.0
        score += 100.0 if promoter_ok else -200.0
        score += 60.0 if hairpin_ok else -120.0

        passed = forbidden_ok and promoter_ok and hairpin_ok
        return score, passed

    def _candidate_score(self, utr: str, codons: list[str], pos: int) -> tuple[float, bool]:
        structural_score, structural_passed = self._local_structural_score(utr, codons, pos)
        codon_score = self._codon_score(codons)
        codons_ok, _, _, _ = self.qualityChecker.codon_checker.run(codons)
        return structural_score + codon_score, structural_passed and codons_ok

    def _better_candidate(
        self,
        trial_score: float,
        trial_passed: bool,
        best_score: float,
        best_passed: bool,
    ) -> bool:
        if trial_passed and not best_passed:
            return True
        if trial_passed == best_passed and trial_score > best_score:
            return True
        return False

    def _optimize_for_rbs(self, peptide: str, utr: str) -> tuple[list[str], dict]:
        """
        Optimize a CDS for one fixed RBS using overlapping codon windows.
        """
        codons = self._seed_codons(peptide) + [self.stop_codon]
        best_codons = codons.copy()
        best_report = self.qualityChecker.run(utr, codons)

        if best_report["passed"]:
            return best_codons, best_report

        windows = self._window_ranges(len(peptide))

        for pass_index in range(self.max_passes):
            changed = False
            ordered_windows = windows if pass_index % 2 == 0 else list(reversed(windows))

            for window_start, window_end in ordered_windows:
                positions = list(range(window_start, window_end))
                if pass_index % 2 == 1:
                    positions.reverse()

                for pos in positions:
                    aa = peptide[pos]
                    current = codons[pos]
                    local_best_codon = current
                    local_best_score, local_best_passed = self._candidate_score(utr, codons, pos)

                    branch_choices = self.codon_choices[aa][: self.max_branch]
                    for alt in branch_choices:
                        if alt == current:
                            continue

                        trial = codons.copy()
                        trial[pos] = alt
                        trial_score, trial_passed = self._candidate_score(utr, trial, pos)

                        if self._better_candidate(
                            trial_score,
                            trial_passed,
                            local_best_score,
                            local_best_passed,
                        ):
                            local_best_codon = alt
                            local_best_score = trial_score
                            local_best_passed = trial_passed

                    if local_best_codon != current:
                        codons[pos] = local_best_codon
                        changed = True

            report = self.qualityChecker.run(utr, codons)
            if self._better_candidate(
                report["score"],
                report["passed"],
                best_report["score"],
                best_report["passed"],
            ):
                best_codons = codons.copy()
                best_report = report

            if best_report["passed"] or not changed:
                break

        return best_codons, best_report

    def run(self, peptide: str, ignores: set) -> Transcript:
        if self.rbsChooser is None:
            self.initiate()

        for aa in peptide:
            if aa not in self.codon_choices:
                raise ValueError(f"Unsupported amino acid: {aa}")

        # Guaranteed fallback so this method never returns None.
        fallback_codons = self._seed_codons(peptide) + [self.stop_codon]
        fallback_cds = "".join(fallback_codons)
        fallback_rbs = self.rbsChooser.run(fallback_cds, ignores)

        best_transcript = Transcript(fallback_rbs, peptide, fallback_codons)
        best_report = self.qualityChecker.run(fallback_rbs.utr, fallback_codons)

        for rbs_option in self._available_rbs_options(ignores):
            try:
                optimized_codons, report = self._optimize_for_rbs(peptide, rbs_option.utr)

                if not optimized_codons:
                    optimized_codons = fallback_codons.copy()
                    report = self.qualityChecker.run(rbs_option.utr, optimized_codons)

                transcript = Transcript(rbs_option, peptide, optimized_codons)

                if self._better_candidate(
                    report["score"],
                    report["passed"],
                    best_report["score"],
                    best_report["passed"],
                ):
                    best_transcript = transcript
                    best_report = report

                if report["passed"]:
                    return transcript

            except Exception:
                continue

        return best_transcript
