from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.transcript_quality_checker import TranscriptQualityChecker


class TranscriptDesigner:
    """
    Reverse-translates a protein sequence into a DNA sequence and chooses
    an RBS while trying to satisfy quality constraints.
    """

    def __init__(self):
        self.codon_choices = {}
        self.rbsChooser = None
        self.qualityChecker = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.qualityChecker = TranscriptQualityChecker()
        self.qualityChecker.initiate()

        # Synonymous codons ordered from stronger to weaker preference.
        self.codon_choices = {
            'A': ["GCG", "GCC", "GCA", "GCT"],
            'C': ["TGC", "TGT"],
            'D': ["GAT", "GAC"],
            'E': ["GAA", "GAG"],
            'F': ["TTC", "TTT"],
            'G': ["GGT", "GGC", "GGA", "GGG"],
            'H': ["CAC", "CAT"],
            'I': ["ATC", "ATT", "ATA"],
            'K': ["AAA", "AAG"],
            'L': ["CTG", "TTA", "TTG", "CTC", "CTA", "CTT"],
            'M': ["ATG"],
            'N': ["AAC", "AAT"],
            'P': ["CCG", "CCA", "CCC", "CCT"],
            'Q': ["CAG", "CAA"],
            'R': ["CGT", "CGC", "AGA", "AGG", "CGA", "CGG"],
            'S': ["TCT", "TCC", "AGC", "TCA", "TCG", "AGT"],
            'T': ["ACC", "ACA", "ACG", "ACT"],
            'V': ["GTG", "GTT", "GTC", "GTA"],
            'W': ["TGG"],
            'Y': ["TAC", "TAT"]
        }

    def _seed_codons(self, peptide: str) -> list[str]:
        return [self.codon_choices[aa][0] for aa in peptide]

    def _choose_best_rbs(self, cds: str, ignores: set):
        return self.rbsChooser.run(cds, ignores)

    def _score(self, utr: str, codons: list[str]) -> float:
        return self.qualityChecker.run(utr, codons)["score"]

    def _optimize_codons(self, peptide: str, utr: str) -> list[str]:
        codons = self._seed_codons(peptide) + ["TAA"]

        # Small greedy search to keep runtime reasonable.
        for _ in range(2):
            improved = False

            for i, aa in enumerate(peptide):
                current = codons[i]
                best_local = current
                best_score = self._score(utr, codons)

                # Try only the top 3 synonyms for speed.
                for alt in self.codon_choices[aa][:3]:
                    if alt == current:
                        continue

                    trial = codons.copy()
                    trial[i] = alt
                    trial_score = self._score(utr, trial)

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
                raise ValueError(f"Unsupported amino acid: {aa}")

        # First-pass CDS for choosing an RBS.
        seed_codons = self._seed_codons(peptide) + ["TAA"]
        seed_cds = "".join(seed_codons)
        selectedRBS = self._choose_best_rbs(seed_cds, ignores)

        # Improve codons using the chosen RBS context.
        optimized_codons = self._optimize_codons(peptide, selectedRBS.utr)

        return Transcript(selectedRBS, peptide, optimized_codons)
