"""Deterministic epoch shuffle + token-budget / study-balanced batching."""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterator, Sequence
from random import Random


def iter_epoch_batches(
    n_items: int,
    *,
    n_tokens: Sequence[int],
    study_ids: Sequence[str],
    task_keys: Sequence[str],
    batch_token_budget: int | None,
    batch_size: int,
    seed: int,
    epoch: int,
) -> Iterator[list[int]]:
    """Yield index lists for one epoch.

    Shuffle is ``Random(seed + epoch)``. Studies are round-robined so a batch
    is not a contiguous pack/study block. Packs until ``batch_token_budget``
    (when set) or ``batch_size`` samples; a single sample may exceed the budget.
    """
    if n_items <= 0:
        return
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    if len(n_tokens) != n_items or len(study_ids) != n_items or len(task_keys) != n_items:
        raise ValueError("n_tokens, study_ids, and task_keys must match n_items")
    rng = Random(int(seed) + int(epoch))  # noqa: S311 — deterministic train shuffle
    order = list(range(n_items))
    rng.shuffle(order)

    by_bucket: dict[str, deque[int]] = defaultdict(deque)
    for idx in order:
        by_bucket[f"{study_ids[idx]}|{task_keys[idx]}"].append(idx)
    bucket_order = list(by_bucket)
    rng.shuffle(bucket_order)
    mixed: list[int] = []
    queues = [by_bucket[s] for s in bucket_order]
    while any(queues):
        for queue in queues:
            if queue:
                mixed.append(queue.popleft())  # noqa: PERF401 — popleft side effect

    budget = int(batch_token_budget) if batch_token_budget else None
    if budget is not None and budget < 1:
        raise ValueError("batch_token_budget must be >= 1 when set")
    max_n = int(batch_size)
    start = 0
    while start < n_items:
        batch: list[int] = []
        tokens = 0
        while start < n_items:
            idx = mixed[start]
            sample_tokens = max(1, int(n_tokens[idx]))
            if batch and (
                len(batch) >= max_n or (budget is not None and tokens + sample_tokens > budget)
            ):
                break
            batch.append(idx)
            tokens += sample_tokens
            start += 1
            if budget is None and len(batch) >= max_n:
                break
        yield batch
