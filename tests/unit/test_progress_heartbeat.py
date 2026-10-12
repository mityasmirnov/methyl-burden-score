"""Heartbeat: quiet when fast, regular when slow, honest ETA, never crashes the loop."""

from __future__ import annotations

from mbs.training.progress import Heartbeat


class _Clock:
    def __init__(self) -> None:
        self.t = 0.0

    def __call__(self) -> float:
        return self.t


def _hb(total: int, interval: float = 300.0) -> tuple[Heartbeat, _Clock, list[str]]:
    clock, lines = _Clock(), []
    return Heartbeat("scoring", total, interval_s=interval, clock=clock, emit=lines.append), clock, lines


def test_fast_loop_stays_silent() -> None:
    hb, clock, lines = _hb(1000)
    for i in range(1, 1001):
        clock.t += 0.01  # 10 s total, well under the 300 s interval
        assert hb.update(i) is False
    hb.finish()
    assert lines == []


def test_slow_loop_emits_once_per_interval_with_correct_eta() -> None:
    hb, clock, lines = _hb(1000, interval=300.0)
    emitted = 0
    for i in range(1, 501):
        clock.t += 1.0  # 1 s/item
        emitted += int(hb.update(i))
    assert emitted == 1  # 500 s elapsed -> only the t=300 tick, not one per call
    assert "300/1000" in lines[0] and "30.0%" in lines[0]
    assert "1.00s/item" in lines[0]
    # 700 items remain at 1 s/item = 700 s ~= 12m
    assert "ETA 12m" in lines[0]


def test_eta_uses_observed_rate_not_a_constant() -> None:
    hb, clock, lines = _hb(100, interval=10.0)
    clock.t = 20.0
    hb.update(10)  # 2 s/item observed; 90 items left -> 180 s -> 3m
    assert "2.00s/item" in lines[0] and "ETA 3m" in lines[0]


def test_zero_done_and_overrun_do_not_crash() -> None:
    hb, clock, lines = _hb(10, interval=1.0)
    clock.t = 5.0
    assert hb.update(0) is False  # would divide by zero
    hb.update(15)  # done > total: remaining clamps at 0
    assert "ETA 0s" in lines[0]


def test_finish_reports_only_long_loops() -> None:
    hb, clock, lines = _hb(10, interval=300.0)
    clock.t = 4000.0  # 66.7 min: under the 90-min hours threshold, so shown in minutes
    hb.finish()
    assert lines and "done 10 in 67m" in lines[-1]
    clock.t = 7200.0
    hb.finish()
    assert "done 10 in 2.0h" in lines[-1]


def test_eta_follows_recent_pace_when_the_job_slows_down() -> None:
    """Regression for the 2026-10-12 scoring pass: ~3x slowdown, ETA stayed optimistic."""
    hb, clock, lines = _hb(1000, interval=300.0)
    done = 0
    for _ in range(300):  # first window: 1 s/item
        clock.t += 1.0
        done += 1
        hb.update(done)
    for _ in range(100):  # second window: 3 s/item
        clock.t += 3.0
        done += 1
        hb.update(done)
    assert len(lines) == 2
    # lifetime average at t=600, done=400 is 1.5 s/item -> would claim 600*1.5 = 15m left;
    # the recent pace is 3 s/item -> 600*3 = 30m left, which is the honest number.
    assert "1.50s/item (recent 3.00)" in lines[1]
    assert "ETA 30m" in lines[1]
