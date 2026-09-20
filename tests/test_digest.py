"""单元测试：圆环片段计算与条带聚类纯算法。"""
from app.digest import circular_fragments, cluster_bands, digest_sample
from app.models import CutPointOut, SampleDigestOut


def _cuts(positions):
    return [CutPointOut(position=p, enzymes=["X"]) for p in positions]


def test_circular_fragments_basic_and_sum_equals_plasmid_length():
    fragments = circular_fragments(_cuts([100, 300, 600]), 1000)
    assert [f.length for f in fragments] == [200, 300, 500]
    assert sum(f.length for f in fragments) == 1000
    assert [f.crosses_origin for f in fragments] == [False, False, True]
    # 跨原点片段：600 -> 100，长度 100 + 1000 - 600
    wrap = fragments[-1]
    assert wrap.start == 600 and wrap.end == 100 and wrap.length == 500


def test_two_cuts_produce_two_fragments():
    fragments = circular_fragments(_cuts([0, 250]), 1000)
    assert [f.length for f in fragments] == [250, 750]
    assert fragments[-1].crosses_origin is True


def test_clustering_gap_larger_than_tolerance_splits_bands():
    # 参考型切点 [100,600] -> 片段 [500,500]；变异型 [110,600] -> [490,510]
    samples = {
        "reference": SampleDigestOut(
            cuts=_cuts([100, 600]),
            fragments=circular_fragments(_cuts([100, 600]), 1000),
            visible_fragment_count=2,
            hidden_fragment_count=0,
        ),
        "variant": SampleDigestOut(
            cuts=_cuts([110, 600]),
            fragments=circular_fragments(_cuts([110, 600]), 1000),
            visible_fragment_count=2,
            hidden_fragment_count=0,
        ),
    }
    # 490/500/500/510 相邻差均 <= 10 -> 单连接聚成 1 条带，两样本各 2 条
    bands = cluster_bands(samples, tolerance=10, min_visible=0)
    assert len(bands) == 1
    assert bands[0].reference_count == 2 and bands[0].variant_count == 2
    assert not bands[0].distinguishable

    # 容差 9 -> 490 | 500,500 | 510：3 条带，全部可区分
    bands_strict = cluster_bands(samples, tolerance=9, min_visible=0)
    assert len(bands_strict) == 3
    assert all(b.distinguishable for b in bands_strict)


def test_min_visible_fragments_are_excluded_from_bands():
    site_map = {
        "A": {"reference": [0, 10], "variant": [0, 10]},
        "B": {"reference": [900], "variant": [900]},
    }
    ref = digest_sample(site_map, "A", "B", "reference", 1000, min_visible=50)
    var = digest_sample(site_map, "A", "B", "variant", 1000, min_visible=50)
    bands = cluster_bands({"reference": ref, "variant": var}, tolerance=20, min_visible=50)
    # 片段为 10 / 890 / 100：10 不可见，890 与 100 相距 790 -> 各成一条带
    assert ref.hidden_fragment_count == 1
    assert len(bands) == 2
