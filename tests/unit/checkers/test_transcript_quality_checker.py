from genedesign.checkers.transcript_quality_checker import TranscriptQualityChecker


def test_transcript_quality_checker_returns_expected_keys():
    checker = TranscriptQualityChecker()
    checker.initiate()

    result = checker.run("AGGAGG", ["ATG", "GCT", "GAA", "TAA"])

    expected_keys = {
        "passed",
        "score",
        "codons_ok",
        "diversity",
        "rare_count",
        "cai",
        "forbidden_ok",
        "forbidden_hit",
        "promoter_ok",
        "promoter_hit",
        "hairpin_ok",
        "hairpin_hit",
    }

    assert expected_keys.issubset(result.keys())
    assert isinstance(result["score"], float)
    assert isinstance(result["passed"], bool)


def test_transcript_quality_checker_detects_forbidden_sequence():
    checker = TranscriptQualityChecker()
    checker.initiate()

    # UTR contains EcoRI site GAATTC, which is forbidden in the existing checker.
    result = checker.run("AGGAGGAATTC", ["ATG", "GCT", "TAA"])

    assert result["forbidden_ok"] is False
    assert result["forbidden_hit"] == "GAATTC"
    assert result["passed"] is False
