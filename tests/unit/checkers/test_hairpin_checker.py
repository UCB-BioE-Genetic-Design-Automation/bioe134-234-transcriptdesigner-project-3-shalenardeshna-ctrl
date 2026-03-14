import genedesign.checkers.hairpin_checker as hc


def test_short_sequence_is_checked(monkeypatch):
    def fake_counter(seq, min_stem, min_loop, max_loop):
        return 2, "BAD_HAIRPIN"

    monkeypatch.setattr(hc, "hairpin_counter", fake_counter)

    ok, hairpin = hc.hairpin_checker("A" * 20)

    assert ok is False
    assert hairpin == "BAD_HAIRPIN"


def test_tail_window_is_checked(monkeypatch):
    def fake_counter(seq, min_stem, min_loop, max_loop):
        if seq.endswith("GGGG"):
            return 2, "TAIL_HAIRPIN"
        return 0, None

    monkeypatch.setattr(hc, "hairpin_counter", fake_counter)

    seq = "A" * 60 + "GGGG"
    ok, hairpin = hc.hairpin_checker(seq)

    assert ok is False
    assert hairpin == "TAIL_HAIRPIN"
