import pandas as pd
import pytest

from prime_fock_nd import (
    NumericBinningSpec,
    CategoricalSpec,
    MultiHotSpec,
    PreEncodedSpec,
    PrimeFockNDConfig,
    PrimeFockNDEncoder,
)


def test_general_encoder_width_and_shape():
    df = pd.DataFrame({
        "time": [1, 1, 2, 2],
        "entity": ["A", "B", "A", "B"],
        "x": [0.1, 2.0, 0.2, 2.2],
        "kind": ["red", "blue", "red", "blue"],
        "tags": ["a|b", "b", "a", "c"],
        "z1": [1.0, 0.0, 0.5, 0.2],
        "z2": [0.0, 1.0, 0.5, 0.0],
    })
    cfg = PrimeFockNDConfig(
        index_cols=["time", "entity"],
        item_col="entity",
        slice_cols=["time"],
        specs=[
            NumericBinningSpec("x", bins=[1.0]),
            CategoricalSpec("kind"),
            MultiHotSpec("tags", separator="|"),
            PreEncodedSpec("latent", columns=["z1", "z2"]),
        ],
    )

    encoded = PrimeFockNDEncoder(cfg).fit_transform(df)
    assert encoded.exponents.shape == (2, 2, 9)
    assert len(encoded.mode_metadata) == 9


def test_duplicate_index_rejected():
    df = pd.DataFrame({
        "time": [1, 1],
        "entity": ["A", "A"],
        "x": [0.1, 0.2],
    })
    cfg = PrimeFockNDConfig(
        index_cols=["time", "entity"],
        item_col="entity",
        slice_cols=["time"],
        specs=[NumericBinningSpec("x", bins=[0.0])],
    )
    enc = PrimeFockNDEncoder(cfg).fit(df)
    with pytest.raises(ValueError, match="Duplicate"):
        enc.transform(df)
