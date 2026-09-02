"""Arm id glossary used in 7G inspection reports."""

from __future__ import annotations

from mbs.inspection.arm_glossary import (
    ARM_DESCRIPTIONS,
    arm_description,
    render_arm_glossary_section,
    render_full_glossary_doc,
)


def test_known_arms_have_descriptions() -> None:
    for arm_id in (
        "P2-G",
        "P4-G",
        "C-mvalue-enet",
        "C-mvalue-enet-G",
        "C-mvalue-enetS",
        "N-cascade-l1",
        "N-cascade-S",
        "N-light-type",
    ):
        desc = arm_description(arm_id)
        assert arm_id in desc or len(desc) > 20
        assert arm_id in ARM_DESCRIPTIONS


def test_render_arm_glossary_section_includes_prefix_legend() -> None:
    md = "\n".join(render_arm_glossary_section(["P2-G", "C-mvalue-enet-G"]))
    assert "## Arm glossary" in md
    assert "Naming prefixes" in md
    assert "`P2-G`" in md
    assert "`C-mvalue-enet-G`" in md


def test_render_full_glossary_doc_covers_eval_modes() -> None:
    doc = render_full_glossary_doc()
    assert "mbs_e2e" in doc
    assert "fusion_full" in doc
    assert "C-mvalue-enetS" in doc
