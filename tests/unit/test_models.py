from __future__ import annotations

import torch

from mbs.models import FlatDeepSet, HierarchicalDeepSet, SeedMaskedLinearHead


def test_flat_model_is_permutation_invariant_and_neutral_for_missing_gene() -> None:
    torch.manual_seed(11)
    model = FlatDeepSet(input_dim=4, dropout=0.0)
    model.eval()

    features = torch.randn(6, 4)
    cpg_to_gene = torch.tensor([0, 0, 1, 1, 1, 0])
    permutation = torch.tensor([5, 2, 0, 4, 1, 3])

    output = model(features, cpg_to_gene, n_gene_instances=3)
    permuted = model(
        features[permutation],
        cpg_to_gene[permutation],
        n_gene_instances=3,
    )

    assert torch.allclose(output["mbs"], permuted["mbs"])
    assert output["present"].tolist() == [True, True, False]
    assert torch.allclose(output["mbs"][2], torch.tensor(0.5))
    assert torch.allclose(output["centered_mbs"][2], torch.tensor(0.0))


def test_hierarchical_model_is_permutation_invariant() -> None:
    torch.manual_seed(17)
    model = HierarchicalDeepSet(
        input_dim=3,
        n_region_types=5,
        dropout=0.0,
    )
    model.eval()

    features = torch.randn(6, 3)
    cpg_to_region = torch.tensor([0, 0, 1, 2, 2, 0])
    region_type = torch.tensor([0, 1, 3, 4])
    region_to_gene = torch.tensor([0, 0, 1, 2])
    permutation = torch.tensor([4, 1, 5, 0, 3, 2])

    output = model(
        cpg_features=features,
        cpg_to_region=cpg_to_region,
        region_type=region_type,
        region_to_gene=region_to_gene,
        n_regions=4,
        n_gene_instances=3,
    )
    permuted = model(
        cpg_features=features[permutation],
        cpg_to_region=cpg_to_region[permutation],
        region_type=region_type,
        region_to_gene=region_to_gene,
        n_regions=4,
        n_gene_instances=3,
    )

    assert torch.allclose(output["mbs"], permuted["mbs"], atol=1e-6)
    assert output["region_present"].tolist() == [True, True, True, False]
    assert output["present"].tolist() == [True, True, False]
    assert torch.allclose(output["mbs"][2], torch.tensor(0.5))


def test_seed_mask_and_presence_control_linear_head() -> None:
    seed_mask = torch.tensor(
        [
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ]
    )
    head = SeedMaskedLinearHead(
        n_genes=3,
        n_outputs=2,
        seed_mask=seed_mask,
    )
    with torch.no_grad():
        head.gene_weight.fill_(2.0)
        head.bias.zero_()

    mbs = torch.tensor([[0.75, 0.25, 0.9]])
    present = torch.tensor([[True, True, False]])
    output = head(mbs, present)

    # Output 0 uses gene 0 and absent gene 2. Output 1 uses gene 1 only.
    assert torch.allclose(output, torch.tensor([[0.5, -0.5]]))
