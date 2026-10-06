# Claims and proof packages

All commands below are run from the workspace produced by `unpack.py`.
The complete per-code routing map is restored as `docs/data-map.json`.

| Claim | Public records | Finite check in expanded workspace |
|---|---|---|
| catalogue frames and upper witnesses | [records](proofs/reference.json) | `python3 verify.py` |
| regular canonical pairing obstruction | [records](proofs/lp_pairing_obstruction.json) | `python3 verification/regular_pairing/check.py` |
| flagship pp width25 | [records](proofs/pp_frame_optimum.json) | `python3 verification/pp_width/check.py --replay-prefix` |
| noncyclic quaternary pp | [records](proofs/noncyclic_pp.json) | `python3 verification/noncyclic_pp/check.py` |
| independent constituent lpq544 | [records](proofs/independent_lp.json) | `python3 verification/lpq_extended/independent544/check.py` |
| lpq1088 distance22 optimal width22 | [records](proofs/lp_replays.json) | `python3 codes/lp-quaternary/lp-f4-1088-128-22-w10-f22/proof/replay.py` |
| sparse check eight lpq | [records](proofs/lp_replays.json) | `python3 codes/lp-quaternary/lp-f4-544-64-14-w8/proof/replay.py` |
| sparse check eight quaternary pp | [records](proofs/check8_pp.json) | `python3 verification/check8_pp/check.py` |
| period48 pp distance20 | [records](proofs/even_period_pp.json) | `python3 verification/even_period_pp/check.py` |
| period48 pp check8 | [records](proofs/check8_pp_period48.json) | `python3 verification/check8_pp_period48/check.py` |
| period48 pp check10 bound22 | [records](proofs/pp48_w10_bound.json) | `python3 verification/pp48_w10_bound/check.py` |
| flagship parent distance24 | [records](proofs/pp_distance.json) | `python3 verification/pp_distance/check.py` |
| five small cyclic pp | [records](proofs/small_quaternary_pp.json) | `python3 verification/small_quaternary_pp/check.py` |

The LP pairing package includes both complete physical orbit pools, pairing
polynomials, compatibility graphs, proper colorings, and separate minimum-weight
sector bases. The PP width package includes the logical detector, attaining
frame, exact engine, and complete exclusion coverage. Large search trees are
not distributed. See [VERIFYING.md](VERIFYING.md) for the distinction between
finite certificates and completed exhaustive executions.
