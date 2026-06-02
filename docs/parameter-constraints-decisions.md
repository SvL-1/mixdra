# Parameter constraints — design decisions

Captured during brainstorming on 2026-05-30. These are the user's answers
verbatim, recording *why* the constraint behaviour is being changed in the
mixdra fitting engine (`R/fit.R`).

## Background

The first binary example dataset returned (CA reference):

| Parameter        | Value       |
|------------------|-------------|
| max              | 0.955598995 |
| beta (Chem 1)    | 0.880510854 |
| beta (Chem 2)    | 2.497718434 |
| EC50 (Chem 1)    | 100         |
| EC50 (Chem 2)    | 0.916824949 |

EC50 (Chem 1) = 100 was a parameter pinned against the **artificial upper
bound** produced by the old seed-relative box constraints
(`upper = start * upper_mult`), not a genuine fit. That motivated rethinking how
parameter constraints work.

## Decisions

### Q1. Is EC50 Chem 1 = 100 a genuine domain bound, or just where our artificial cap landed?

> The model should sometimes be bound to biological or chemical constraints
> beyond which extrapolation makes no sense, like a concentration above 100% or
> constraining the beta (slope) of the curve to a certain maximum, as biological
> relevance of extreme curve slopes goes down, while requiring substantial effort
> to optimize. Do you understand what I mean?

**Takeaway:** bounds should encode domain meaning (a concentration cannot exceed
its physically meaningful range; extreme slopes are biologically implausible and
costly to optimise), not be inferred from the seed value.

### Q2. How should the default bounds be determined when the user doesn't specify them explicitly?

> All constraints should be implemented manually by the user, except for all
> params get only positivity. Also, a and b should remain unconstrained, as these
> values are used in later analysis of e.g. the concentration at which the
> interactions switch from synergistic to antagonistic.

**Takeaway:**
- Constraints are **manual** — the user supplies them per parameter.
- Default for base curve params (`max`, `slope*`, `ec50*`): **positivity only**
  (lower `~1e-8`, upper `Inf`). The old seed-multiplier scheme is removed.
- Deviation parameters `a` / `b` stay **unconstrained**, because they feed
  downstream interaction analysis (e.g. the concentration at which the
  interaction switches from synergistic to antagonistic). Bounding them would
  distort that analysis. (The ternary per-chemical `b1/b2/b3` of dose-ratio were
  removed on 2026-06-02 when DR/DL were restricted to binary mixtures only.)

### Correctness exception (raised in design, pending confirmation)

For **binary** data, `max` is the control-response *probability*; the model
yields fitted probabilities > 1 if `max` exceeds 1. Binary `max` therefore keeps
a hard default `upper = 1` for correctness (the user may still tighten it or fix
it, e.g. Sam fixes it to 0.98).
