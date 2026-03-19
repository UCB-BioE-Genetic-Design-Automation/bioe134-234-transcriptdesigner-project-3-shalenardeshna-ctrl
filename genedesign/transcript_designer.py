from collections import defaultdict

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.transcript_quality_checker import TranscriptQualityChecker


class TranscriptDesigner:
    """
    Reverse-translates a protein sequence into a DNA sequence and chooses
    an RBS while trying to satisfy quality constraints.

    This version keeps the original model interfaces, but adds a more focused
    hairpin-repair stage so benchmark hairpin failures are actively targeted.
    """

    def __init__(self):
        self.codon_choices = {}
        self.rbsChooser = None
        self.qualityChecker = None

        # Balanced for speed + stronger hairpin cleanup.
        self.window_codons = 12
        self.window_step = 6
        self.max_passes = 2
        self.max_seed_phases = 2
        self.max_rbs_options = 2
        self.max_hairpin_cleanup_rounds = 2
        self.max_hairpin_positions = 12
        self.stop_codon = "TAA"

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.qualityChecker = TranscriptQualityChecker()
        self.qualityChecker.initiate()

        self._init_codon_choices()

    def _init_codon_choices(self) -> None:
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

        freqs = self.qualityChecker.codon_checker.codon_frequencies

        self.codon_choices = {}
        for aa, codons in genetic_code.items():
            ordered = codons[:]
            ordered.sort(key=lambda codon: freqs.get(codon, 0.0), reverse=True)
            self.codon_choices[aa] = ordered

    def _available_rbs_options(self, ignores: set):
        options = [opt for opt in self.rbsChooser.rbsOptions if opt not in ignores]
        if not options:
            raise Exception("No valid RBS option available.")
        return options

    def _better_report(self, new_report: dict, old_report: dict | None) -> bool:
        if old_report is None:
            return True
        if new_report["rank"] < old_report["rank"]:
            return True
        if new_report["rank"] == old_report["rank"] and new_report["score"] > old_report["score"]:
            return True
        return False

    def _seed_codons(self, peptide: str, phase: int = 0) -> list[str]:
        counts = defaultdict(int)
        codons = []

        for i, aa in enumerate(peptide):
            choices = self.codon_choices[aa]

            if i == 0 and aa == "M":
                codons.append("ATG")
                counts[aa] += 1
                continue

            idx = (counts[aa] + phase) % len(choices)
            codons.append(choices[idx])
            counts[aa] += 1

        return codons

    def _window_ranges(self, n_positions: int) -> list[tuple[int, int]]:
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

    def _nt_range_to_codon_positions(
        self,
        nt_start: int,
        nt_end: int,
        utr_len: int,
        peptide_len: int,
    ) -> list[int]:
        cds_start = max(0, nt_start - utr_len)
        cds_end = max(0, nt_end - utr_len)

        first = cds_start // 3
        last_exclusive = min(peptide_len, (cds_end + 2) // 3)

        if last_exclusive <= first:
            return []

        return list(range(first, last_exclusive))

    def _add_motif_positions(
        self,
        positions: set[int],
        transcript_dna: str,
        motif: str | None,
        utr_len: int,
        peptide_len: int,
    ) -> None:
        if not motif:
            return

        motif = str(motif).upper()
        if not motif:
            return

        start = transcript_dna.find(motif)
        while start != -1:
            end = start + len(motif)
            for pos in self._nt_range_to_codon_positions(start, end, utr_len, peptide_len):
                positions.add(pos)
            start = transcript_dna.find(motif, start + 1)

    def _hairpin_pressure(
        self,
        peptide: str,
        utr: str,
        report: dict,
    ) -> dict[int, int]:
        pressure = {}
        utr_len = len(utr)
        peptide_len = len(peptide)

        for detail in report.get("hairpin_bad_window_details", []):
            weight = max(1, detail.get("count", 2) - 1)
            for pos in self._nt_range_to_codon_positions(
                detail["start"],
                detail["end"],
                utr_len,
                peptide_len,
            ):
                pressure[pos] = pressure.get(pos, 0) + weight

        return pressure

    def _problem_positions(
        self,
        peptide: str,
        utr: str,
        codons: list[str],
        report: dict,
    ) -> list[int]:
        positions = set()
        transcript_dna = (utr + "".join(codons)).upper()
        pressure = self._hairpin_pressure(peptide, utr, report)

        for pos in pressure:
            positions.add(pos)

        self._add_motif_positions(
            positions,
            transcript_dna,
            report.get("forbidden_hit"),
            len(utr),
            len(peptide),
        )
        self._add_motif_positions(
            positions,
            transcript_dna,
            report.get("promoter_hit"),
            len(utr),
            len(peptide),
        )

        if not report.get("codons_ok", True) or not positions:
            for pos, aa in enumerate(peptide):
                if len(self.codon_choices[aa]) > 1:
                    positions.add(pos)

        ordered_positions = sorted(
            positions,
            key=lambda pos: (-pressure.get(pos, 0), pos),
        )
        return ordered_positions

    def _ordered_candidate_codons(self, aa: str, current: str, phase: int) -> list[str]:
        choices = self.codon_choices[aa][:]
        if len(choices) <= 1:
            return []

        shift = phase % len(choices)
        rotated = choices[shift:] + choices[:shift]

        return [codon for codon in rotated if codon != current]

    def _repair_once(
        self,
        peptide: str,
        utr: str,
        codons: list[str],
        report: dict,
        phase: int,
    ) -> tuple[list[str], dict]:
        current_codons = codons[:]
        current_report = report
        target_positions = self._problem_positions(peptide, utr, current_codons, current_report)
        windows = self._window_ranges(len(peptide))

        prioritized_windows = []
        for start, end in windows:
            if any(start <= pos < end for pos in target_positions):
                prioritized_windows.append((start, end))
        for window in windows:
            if window not in prioritized_windows:
                prioritized_windows.append(window)

        for window_index, (start, end) in enumerate(prioritized_windows):
            window_positions = [pos for pos in target_positions if start <= pos < end]
            if not window_positions:
                continue

            for pos in window_positions:
                aa = peptide[pos]
                current = current_codons[pos]

                best_local_codons = current_codons
                best_local_report = current_report

                candidate_phase = phase + window_index + pos
                for alt in self._ordered_candidate_codons(aa, current, candidate_phase):
                    trial_codons = current_codons[:]
                    trial_codons[pos] = alt
                    trial_report = self.qualityChecker.run(utr, trial_codons)

                    if self._better_report(trial_report, best_local_report):
                        best_local_codons = trial_codons
                        best_local_report = trial_report

                if self._better_report(best_local_report, current_report):
                    current_codons = best_local_codons
                    current_report = best_local_report

                    if current_report["passed"]:
                        return current_codons, current_report

        return current_codons, current_report

    def _focused_hairpin_cleanup(
        self,
        peptide: str,
        utr: str,
        codons: list[str],
        report: dict,
    ) -> tuple[list[str], dict]:
        current_codons = codons[:]
        current_report = report

        for cleanup_round in range(self.max_hairpin_cleanup_rounds):
            if current_report.get("hairpin_total_excess", 0) <= 0:
                break

            pressure = self._hairpin_pressure(peptide, utr, current_report)
            if not pressure:
                break

            candidate_positions = [
                pos
                for pos, _ in sorted(pressure.items(), key=lambda item: (-item[1], item[0]))
            ][: self.max_hairpin_positions]

            improved = False

            for pos in candidate_positions:
                aa = peptide[pos]
                current = current_codons[pos]

                best_local_codons = current_codons
                best_local_report = current_report

                for alt in self._ordered_candidate_codons(aa, current, cleanup_round + pos):
                    trial_codons = current_codons[:]
                    trial_codons[pos] = alt
                    trial_report = self.qualityChecker.run(utr, trial_codons)

                    if self._better_report(trial_report, best_local_report):
                        best_local_codons = trial_codons
                        best_local_report = trial_report

                if self._better_report(best_local_report, current_report):
                    current_codons = best_local_codons
                    current_report = best_local_report
                    improved = True

                    if current_report["passed"]:
                        return current_codons, current_report

            if not improved:
                break

        return current_codons, current_report

    def _optimize_for_rbs(self, peptide: str, rbs_option) -> tuple[list[str], dict]:
        best_codons = None
        best_report = None

        max_phase = max(
            1,
            min(
                self.max_seed_phases,
                max(len(choices) for choices in self.codon_choices.values()),
            ),
        )

        for phase in range(max_phase):
            codons = self._seed_codons(peptide, phase=phase) + [self.stop_codon]
            report = self.qualityChecker.run(rbs_option.utr, codons)

            if self._better_report(report, best_report):
                best_codons = codons[:]
                best_report = report

            if report["passed"]:
                return codons, report

            working_codons = codons[:]
            working_report = report

            for repair_round in range(self.max_passes):
                new_codons, new_report = self._repair_once(
                    peptide,
                    rbs_option.utr,
                    working_codons,
                    working_report,
                    phase + repair_round,
                )

                if self._better_report(new_report, best_report):
                    best_codons = new_codons[:]
                    best_report = new_report

                if new_report["passed"]:
                    return new_codons, new_report

                if (
                    new_report["rank"] == working_report["rank"]
                    and new_report["score"] <= working_report["score"]
                ):
                    break

                working_codons = new_codons
                working_report = new_report

            if best_report is not None and not best_report["hairpin_ok"]:
                cleaned_codons, cleaned_report = self._focused_hairpin_cleanup(
                    peptide,
                    rbs_option.utr,
                    best_codons,
                    best_report,
                )
                if self._better_report(cleaned_report, best_report):
                    best_codons = cleaned_codons[:]
                    best_report = cleaned_report

                if cleaned_report["passed"]:
                    return cleaned_codons, cleaned_report

        return best_codons, best_report

    def run(self, peptide: str, ignores: set) -> Transcript:
        if self.rbsChooser is None:
            self.initiate()

        for aa in peptide:
            if aa not in self.codon_choices:
                raise ValueError(f"Unsupported amino acid: {aa}")

        fallback_codons = self._seed_codons(peptide, phase=0) + [self.stop_codon]
        fallback_rbs = self.rbsChooser.run("".join(fallback_codons), ignores)
        best_transcript = Transcript(fallback_rbs, peptide, fallback_codons)
        best_report = self.qualityChecker.run(fallback_rbs.utr, fallback_codons)

        for rbs_option in self._available_rbs_options(ignores)[: self.max_rbs_options]:
            try:
                codons, report = self._optimize_for_rbs(peptide, rbs_option)
                transcript = Transcript(rbs_option, peptide, codons)

                if self._better_report(report, best_report):
                    best_transcript = transcript
                    best_report = report

                if report["passed"]:
                    return transcript
            except Exception:
                continue

        return best_transcript
