"""Platform id aliases (Hub sample-info → catalog / matrix manifests)."""

from __future__ import annotations

import pandas as pd

PLATFORM_ALIASES = {
    "450K": "HM450",
    "450k": "HM450",
    "HM450": "HM450",
    "EPIC": "EPIC",
    "EPICv2": "EPICv2",
}


def normalize_platform(raw: object | None) -> str | None:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    text = str(raw).strip()
    if text == "":
        return None
    return PLATFORM_ALIASES.get(text, text)
