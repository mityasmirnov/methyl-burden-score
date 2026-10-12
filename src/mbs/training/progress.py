"""Time-based heartbeat for long per-sample loops.

On 2026-10-07 a full-width run sat in post-training scoring for 4h24m with no log
output, holding ~43 GiB of GPU memory at 0% utilisation, and nobody could tell
whether it was minutes from done or hours. Per-sample passes (orientation, test
pass, Stage A scoring) have no natural epoch boundary to log at, so this emits at
most one line per ``interval_s`` seconds with a running ETA. Time-based rather than
count-based so fast jobs (unit tests, 65k smokes) stay quiet and slow ones stay
visible.
"""

from __future__ import annotations

import time
from collections.abc import Callable


def _fmt_duration(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    if seconds < 90:
        return f"{seconds:.0f}s"
    if seconds < 5400:
        return f"{seconds / 60:.0f}m"
    return f"{seconds / 3600:.1f}h"


class Heartbeat:
    """Emit ``[progress] <label> done/total (pct) rate (recent) ETA`` at most every ``interval_s``.

    ETA is computed from the rate over the last window, not the lifetime average."""

    def __init__(
        self,
        label: str,
        total: int,
        *,
        interval_s: float = 300.0,
        clock: Callable[[], float] = time.monotonic,
        emit: Callable[[str], None] | None = None,
    ) -> None:
        self.label = label
        self.total = int(total)
        self.interval_s = float(interval_s)
        self._clock = clock
        self._emit = emit if emit is not None else self._print
        self._t0 = clock()
        self._last = self._t0
        self._last_done = 0

    @staticmethod
    def _print(msg: str) -> None:
        print(msg, flush=True)  # noqa: T201

    def update(self, done: int) -> bool:
        """Report progress if the interval has elapsed. Returns True if a line was emitted."""
        now = self._clock()
        if now - self._last < self.interval_s or done <= 0:
            return False
        elapsed = now - self._t0
        per_item = elapsed / done
        # ETA must follow the CURRENT pace, not the lifetime average: on 2026-10-12 a
        # scoring pass slowed ~3x partway through and a cumulative-average ETA kept
        # printing "24m" while the real remaining time was roughly double. Rate is
        # measured over the window since the previous line.
        window_items = done - self._last_done
        recent = (now - self._last) / window_items if window_items > 0 else per_item
        self._last, self._last_done = now, done
        remaining = max(self.total - done, 0) * recent
        pct = 100.0 * done / self.total if self.total else 100.0
        self._emit(
            f"[progress] {self.label} {done}/{self.total} ({pct:.1f}%) "
            f"{per_item:.2f}s/item (recent {recent:.2f}) elapsed {_fmt_duration(elapsed)} "
            f"ETA {_fmt_duration(remaining)}"
        )
        return True

    def finish(self) -> None:
        """Emit one closing line with the total wall time (only for loops that ran long)."""
        elapsed = self._clock() - self._t0
        if elapsed >= self.interval_s:
            self._emit(f"[progress] {self.label} done {self.total} in {_fmt_duration(elapsed)}")
