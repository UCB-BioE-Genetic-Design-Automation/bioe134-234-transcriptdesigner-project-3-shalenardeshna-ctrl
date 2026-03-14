from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker


class TranscriptDesigner:
    """
    Reverse-translates a peptide into a CDS while trying to satisfy:
    - high CAI
    - low hairpin count
    - no forbidden sequences
    - no internal promoters

    Keeps the original Transcript model unchanged.
    """

    def __init__(self):
        self.codon_choices = {}
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenChecker = None
        self.promoterChecker = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        # Synonymous codons ordered roughly from stronger to weaker E. coli preference.
        self.codon_choices = {
            "A": ["GCG", "GCC", "GCA", "GCT"],
            "C": ["TGC", "TGT"],
            "D": ["GAT", "GAC"],
            "E": ["GAA", "GAG"],
            "F": ["TTC", "TTT"],
            "G": ["GGT", "GGC", "GGA", "GGG"],
            "H": ["CAC", "CAT"],
            "I": ["ATC", "ATT", "ATA"],
            "K": ["AAA", "AAG"],
            "L": ["CTG", "TTA", "TTG", "CTC", "CTA", "CTT"],
            "M": ["ATG"],
            "N": ["AAC", "AAT"],
            "P": ["CCG", "CCA", "CCC", "CCT"],
            "Q": ["CAG", "CAA"],
            "R": ["CGT", "CGC", "AGA", "AGG", "CGA", "CGG"],
            "S": ["TCT", "TCC", "AGC", "TCA", "TCG", "AGT"],
            "T": ["ACC", "ACA", "ACG", "ACT"],
            "V": ["GTT", "GTG", "GTC", "GTA"],
            "W": ["TGG"],
            "Y": ["TAC", "TAT"],
        }

    def _score_candidate(self, codons_no_stop: list[str]) -> float:
        """
        Higher score is better.
        """
        cds = "".join(codons_no_stop) + "TAA"

        codon_ok, diversity, rare_count, cai = self.codonChecker.run(codons_no_stop)
        forbidden_ok, _ = self.forbiddenChecker.run(cds)
        promoter_ok, _ = self.promoterChecker.run(cds)
        hairpin_ok, _ = hairpin_checker(cds)

        score = 0.0

        # Hard constraints get large rewards/penalties.
        score += 100.0 if forbidden_ok else -200.0
        score += 100.0 if promoter_ok else -200.0
        score += 60.0 if hairpin_ok else -120.0
        score += 20.0 if codon_ok else -20.0

        # Soft preferences.
        score += 40.0 * cai
        score += 10.0 * diversity
        score -= 8.0 * rare_count

        return score

    def _seed_codons(self, peptide: str) -> list[str]:
        """
        Start from the highest-preference codon for each amino acid.
        """
        return [self.codon_choices[aa][0] for aa in peptide]

    def _optimize_codons(self, peptide: str) -> list[str]:
        """
        Greedy local search over synonymous codons.
        """
        codons = self._seed_codons(peptide)

        for _ in range(4):  # a few passes is enough for this assignment
            improved = False

            for i, aa in enumerate(peptide):
                current = codons[i]
                best_local = current
                best_score = self._score_candidate(codons)

                for alt in self.codon_choices[aa]:
                    if alt == current:
                        continue

                    trial = codons.copy()
                    trial[i] = alt
                    trial_score = self._score_candidate(trial)

                    if trial_score > best_score:
                        best_score = trial_score
                        best_local = alt

                if best_local != current:
                    codons[i] = best_local
                    improved = True

            if not improved:
                break

        return codons

    def run(self, peptide: str, ignores: set) -> Transcript:
        if self.rbsChooser is None:
            self.initiate()

        for aa in peptide:
            if aa not in self.codon_choices:
                raise ValueError(f"Unsupported amino acid in peptide: {aa}")

        codons_no_stop = self._optimize_codons(peptide)
        codons = codons_no_stop + ["TAA"]
        cds = "".join(codons)

        selectedRBS = self.rbsChooser.run(cds, ignores)

        return Transcript(selectedRBS, peptide, codons)


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    print(transcript)
