# Selective correction — Stage A implementation record

Date: 2026-09-09 (Asia/Shanghai)

This record covers only the authorized Stage A CPU arithmetic implementation
and the required legacy regressions. No Habitat/runtime/data preparation,
download, install, clone, training, GPU work, Stage B, pressure grid, or
continuous synthetic band was run. The contract and production operator were
not modified.

## Implementation boundary

The same `idea-stage/conductance_star_arithmetic_check.js` path now has an
explicit `--selective-correction-stage-a` entry. No flag still selects the old
toy-training path. The Stage A path reuses the existing `starFlow`, implicit
cell/Newton solve, arrowhead response, three route credits, normalized
readout, pending/consume semantics, and fixed four-parameter/five-state
fixture.

The new conditional path implements, for fixed incoming during the whole
implicit solve,

```text
q_j = (z_j_minus - s_minus)^2
P   = mean(q)
e_j = (q_j - P) / Cstar
p   = (1, -1, -1, 1)
w_j = 1 + p_j * e_j
alpha_j = k * exp(theta_j * w_j)
I_j = k * (z_j - s) + alpha_j * (z_j - s)^3
```

The derivative carries the same frozen `w` in all three route credits. The
detached diagnostic uses the same forward formula and only removes the soma
error response in the credit path. The common-only mode is the preserved
`w=1` original-R coefficient path. Cold candidates retain `w=1`, and their A
and B route derivatives are retained rather than gated to zero. The code
records `q/P/Cstar/e/w/alpha/differentialConductance`, all three route
credits, full/detached gradients, finite differences, and each theta-version
event without storing the old 12.8 MB structure report inside the new JSON.

## Source and artifact hashes

| Item | SHA-256 |
|---|---|
| JS before Stage A, contract reference | `A6A65444EE7BBF4F08958283C3655FB2931AAF2043BAE67AEEB211F905FC5786` |
| JS after implementation | `FB3B7088A7A286FE99DFD8A68EA245851FE20FF67E3CA0AA38741398724D7629` |
| `ttt02_quick/operator.py` (unchanged) | `C83F26B3233088B72098C5156CCA88DF0ABA0ED4FFB2CC2373F93F43777579A8` |
| Stage A JSON | `A874F4CB7C4F6FD51A72BC4E21E4B9777CD9E1D9C525E8670259493C2A278CBC` |
| preserved constructive JSON | `4A0AA4CAFD2310E0DCD6EA66FB718BD779CA05F0B8596A4A964B069C7A7609EB` |
| preserved association JSON | `1CFC8B4F47261F46F4CF16BFDF98A69301D502A0F2705273249030375192C2A4` |
| preserved S JSON | `CE124DC26B855A5E96EB91BB0F812D89BC3C5BDA8839C4CF4B2F8C20691E6E56` |

The Stage A JSON is 2,271,531 bytes. Its final stdout was also 2,271,531
UTF-8 bytes with the same SHA-256.

## Commands and actual execution results

The Node executable for every command was:

`C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe`

| Command | Exit | Raw result |
|---|---:|---|
| `node.exe --check D:/EV-TTT/idea-stage/conductance_star_arithmetic_check.js` | 0 | empty stdout/stderr; empty-stream SHA `E3B0C44298FC1C149AFBF4C8996FB92427AE41E4649B934CA495991B7852B855` |
| `node.exe D:/EV-TTT/idea-stage/conductance_star_arithmetic_check.js --selective-correction-stage-a` | 2 | JSON written; stdout byte-identical to JSON; exit 2 is the frozen scientific-gate failure, not a runtime crash |
| `node.exe D:/EV-TTT/idea-stage/conductance_star_arithmetic_check.js --constructive-coupling-only` | 0 | 78,490 bytes; stdout/file byte-identical; 132 checks |
| `node.exe D:/EV-TTT/idea-stage/conductance_star_arithmetic_check.js --association-credit-kill-only` | 0 | 67,772 bytes; stdout/file byte-identical; 118 checks |
| `node.exe D:/EV-TTT/idea-stage/conductance_star_arithmetic_check.js --structure-computation-only` | 0 | 12,787,777 bytes; stdout/file byte-identical; preserved S SHA above |

The three legacy files were not replaced by the Stage A command. The explicit
structure-only regression rewrote its own same-path artifact and reproduced
the saved bytes exactly.

## Failed implementation attempts retained in the execution record

These failures occurred before the final Stage A JSON and were fixed without
changing a scientific threshold or adding a fallback:

1. The first post-edit syntax check exited 1 at the output ternary because
   `: output` was missing. It was corrected; the final syntax check is exit 0.
2. The first Stage A launch exited 1 with
   `ReferenceError: stageFiniteDifferenceSteps is not defined`; the declared
   name was `stageAFiniteDifferenceSteps`.
3. The next launch exited 1 when a non-finite coefficient reached `starFlow`.
   The explicit diagnostic showed the update state was already non-finite.
4. The next launch showed the cause was the equal-increment zero-evidence
   check reading an undefined pending `margin`, producing `pBefore/pPlus/delta`
   as non-finite. The Stage A pending cache was corrected to retain
   `margin: computed.value`, as required by the existing evidence interface.
5. The following launch exited 1 on an output-name typo,
   `ReferenceError: deltaLogAlpha is not defined`; it was corrected to use the
   computed `logAlphaDelta` value.

No failed attempt produced a result JSON. The final exit-2 result is the
authoritative raw Stage A artifact.

## Stage A raw gate result

The final JSON records 3 basepoints, 4 evidence relations, 2 orders, 24
two-event units, and both `eta=.02` and `rho=.001` budgets per unit. It has
four reported modes: `full`, `differenceOnly`, `detachedError`, and
`commonOnlyOriginalR`.

Correctness and mechanism outcomes are kept separate:

```text
correctness checks: 10/10 passed
scientific checks:  22/27 passed
basepoint gates:     2/3 passed
24-unit gates:      20/24 passed
Cstar expected:      0.019658382399647234
Cstar recomputed:    0.019658382399647234
stageB implemented:  false
stageB eligible:     false
```

Basepoint full-panel values were:

| Basepoint | `S_cross` | `c` | minimum normalized-Gram eigenvalue | gate |
|---|---:|---:|---:|---|
| cold candidates | `0.01978055397389469` | `0.019331602225113233` | `0.9806683977748868` | pass |
| warm candidates H0/H1 | `0.5397917579915925` | `-0.3736793288992411` | `0.626320671100759` | fail |
| warm candidates H1/H0 | `0.14316644510282872` | `0.10629989516590697` | `0.893700104834093` | pass |

The four failed units are exactly the `y+`,`y+` and `y-`,`y-` same-sign
relations at the warm H0/H1 candidate assignment, in both consumption orders.
Their `U` and `D` values remain positive, and the rho first-B comparisons are
above the original-R 25% requirement (ratios approximately 1.985–2.000), but
the frozen interference/retention thresholds fail. Representative worst
raw values are:

```text
eta H0 same-sign: H/U = 0.5362884334362544, D/U = 0.46371156656374557
eta H1 same-sign: H/U = 0.26037528121433257, D/U = 0.7396247187856674
rho H0 same-sign: H/U = 0.37404990417869644, D/U = 0.6259500958213036
rho H1 same-sign: H/U = 0.373659009163818,   D/U = 0.626340990836182
```

The `H/U <= .25` and `D/U >= .50` failures are reported as failures; no
threshold, fixture, order, relation, or budget was altered. The JSON retains
all four failed units and all passing units, including their two events,
three-route gradients, fixed incoming weight data, actual readout changes,
first-order readout-dot-`DeltaTheta` predictions, and rho reference values.

## Disposition

Stage A arithmetic is implemented and fully reported, but the first frozen
scientific gate did not pass. Stage B is not implemented or started. The
correct next action is Astra analysis of the warm H0/H1 same-sign and
candidate-incoming failure; no formula redesign, threshold relaxation,
pressure expansion, or visual claim is made here.
