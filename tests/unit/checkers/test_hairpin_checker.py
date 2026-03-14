import genedesign.checkers.hairpin_checker as hairpin_checker_module


def test_hairpin_checker_returns_true_when_no_bad_hairpins(monkeypatch):
    calls = []

    def fake_hairpin_counter(chunk, min_stem, min_loop, max_loop):
        calls.append((chunk, min_stem, min_loop, max_loop))
        return 0, None

    monkeypatch.setattr(hairpin_checker_module, "hairpin_counter", fake_hairpin_counter)

    dna = "A" * 100
    passed, hairpin = hairpin_checker_module.hairpin_checker(dna)

    assert passed is True
    assert hairpin is None
    assert len(calls) == 3  # windows starting at 0, 25, 50
    assert all(len(call[0]) == 50 for call in calls)
    assert all(call[1:] == (3, 4, 9) for call in calls)


def test_hairpin_checker_returns_false_when_chunk_has_more_than_one_hairpin(monkeypatch):
    calls = []

    def fake_hairpin_counter(chunk, min_stem, min_loop, max_loop):
        calls.append(chunk)
        if len(calls) == 2:
            return 2, "bad_hairpin"
        return 0, None

    monkeypatch.setattr(hairpin_checker_module, "hairpin_counter", fake_hairpin_counter)

    dna = "A" * 100
    passed, hairpin = hairpin_checker_module.hairpin_checker(dna)

    assert passed is False
    assert hairpin == "bad_hairpin"
    assert len(calls) == 2  # should stop as soon as it finds a bad chunk


def test_hairpin_checker_uses_50bp_windows_with_25bp_overlap(monkeypatch):
    observed_chunks = []

    def fake_hairpin_counter(chunk, min_stem, min_loop, max_loop):
        observed_chunks.append(chunk)
        return 0, None

    monkeypatch.setattr(hairpin_checker_module, "hairpin_counter", fake_hairpin_counter)

    dna = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "ABCDEFGHIJKLMNOPQRSTUV"
    )  # length 100

    passed, hairpin = hairpin_checker_module.hairpin_checker(dna)

    assert passed is True
    assert hairpin is None

    expected_chunks = [
        dna[0:50],
        dna[25:75],
        dna[50:100],
    ]

    assert observed_chunks == expected_chunks
