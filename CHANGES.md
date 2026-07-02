# CHI Heuristic Improvements — Summary of Changes

This document summarizes all algorithmic changes made to the CHI heuristic
during the optimization process, for reference in the paper.

## 1. Pre-Replication of Bestselling SKUs (`REPLICATEALL!`)

**File:** `functions/functions_chisquare.jl` (new function, two methods: sparse + dense Q)
**Called from:** `heuristics/heuristic_chisquare.jl`, before the Phase 1 allocation loop

**What:** Before the dependency-aware Phase 1 allocation loop, the most frequently
ordered SKUs are replicated to ALL D warehouses. This consumes the "bestselling
budget" — the extra capacity beyond what is needed to store each SKU once.

**Budget formula:** `N_tilde / (D-1)` where `N_tilde = sum(capacity) - sum(sku_weight)`,
computed via the existing `BESTSELLING_SKUS()` function from Catalan & Fisher (2012).

**Ranking:** SKUs are ranked by `sum_nor` (the independent coappearance sum
t^E, per unit weight), with transaction frequency as tie-break — the same
order used for the Phase 1 allocation (`nor_order`). An earlier intermediate
version ranked the replication by pure transaction frequency instead, arguing
that `sum_nor` inflates SKUs appearing in large orders. A paired A/B test
(2026-06-11; 2,000 SKUs, 20,000 orders, all six dependency configs, K in
{2, 4, 10}, 20% buffer, 5 seeds) showed the t^E ranking is the better choice:
it wins clearly on uniform-frequency dependency data (up to ~7 pp lower split
ratio on HD-SF, ~2-4 pp on MD-SF) and is statistically tied everywhere else
(frequency's only edge is ~0.1 pp on ID-VF). Intuition: the replication
budget should prioritize SKUs whose coappearance mass cannot be captured by
co-locating a dependency cluster; t^E excludes the dependency premium and
does exactly that.

**Weight-awareness:** Each SKU's `sku_weight[i]` is checked against remaining
warehouse capacity before placement. The budget is consumed in weight units.
This extends EMCI's replication approach (which requires uniform weights) to
heterogeneous SKU weights.

**Degeneration:** When `N_tilde <= 0` (no extra capacity, i.e., buffer = 0), the
function returns immediately. The algorithm is identical to the original CHI
heuristic in this case.

**State updates:** After replication, the function updates `X`, `cap_left`,
`allocated`, `n_allocated`, `sum_dep`, `sum_nor`, `state_dep`, and `state_nor`
so that Phase 1 sees the correct post-replication state.

## 2. Warehouse Weight Recomputation After Replication

**File:** `heuristics/heuristic_chisquare.jl`, after `REPLICATEALL!` call

**What:** After pre-replication changes the capacity landscape, `WHWEIGHT()` is
called again with the updated `cap_left` (rounded to Int64) and the current
`sum_nor` (where replicated SKUs have been zeroed). This ensures that `SELECTIK`
in Phase 1 uses accurate warehouse density weights that reflect the
post-replication capacity distribution.

**Why:** Without this recomputation, SELECTIK's warehouse selection uses weights
computed for the original (pre-replication) capacity, which can misguide
allocation decisions — particularly when replication consumed significant
capacity from all warehouses uniformly.

## 3. Relative Effect Size Filter in Hypothesis Test

**File:** `functions/functions_chisquare.jl`, in `HYOPTHESISTEST_SPARSE()`
**Parameter:** `min_effect::Float64` (default 0.0), passed through from `CHISQUAREHEUR`

**What:** After the chi-square test determines statistical significance, an
additional practical significance check filters out dependencies where the
excess coappearance is less than `min_effect` fraction of the expected
independent coappearance:

```julia
if chi > accept
    val = qij - independent
    if val >= independent * min_effect
        # record dependency
    end
end
```

**Why:** The Bonferroni-corrected chi-square test controls statistical
significance but not practical significance. On variable-frequency (Zipf) data,
popular SKU pairs have large expected coappearances. Random fluctuations of 1-2
extra coappearances can pass the chi-square threshold while representing a tiny
fraction (< 0.3%) of the expected value. These false positives cause
`ADDDEPENDENT!` to pull popular SKUs into the same warehouse unnecessarily.

**Scaling property:** The filter naturally adapts to the frequency structure:
- Popular SKU pairs (high `independent`): need larger absolute excess to pass
- Rare SKU pairs (low `independent`): even small excess passes easily
- Real dependencies (any strength): excess is typically 5-275% of expected
- False positives on popular pairs: excess is typically < 1% of expected

**Recommended value:** `min_effect = 0.01` (1%). This filters noise on popular
pairs without damaging real dependencies on high/medium-dependency data.

## 4. Optimized Significance Level

**File:** `main_benchmark_settings.jl`

**What:** Changed `sig_levels` from `[1.0e-2]` to `[1.0e-5]`.

**Why:** Systematic comparison across six alpha levels (0.01, 0.001, 0.0001,
1e-5, 1e-7, 1e-8) showed that alpha = 1e-5 minimizes average parcels across
all dependency types. It provides the best balance between detecting real
dependencies (HD/MD scenarios) and suppressing false positives (VF scenarios).
The Bonferroni correction dominates the effective threshold regardless of alpha,
but the stricter base level provides additional false-positive control.

## 5. Type Fix in REMOVEALLOC

**File:** `functions/functions_chisquare.jl`, in `REMOVEALLOC()`

**What:** Changed `dep_new = copy(Q)` to
`dep_new = convert(SparseMatrixCSC{Float64,Int64}, copy(Q))`.

**Why:** When `Q` is `SparseMatrixCSC{Int64}` and `dep` contains `Float64`
values, writing `qval - current` (Float64) into the Int64 sparse matrix caused
`InexactError`. This only manifested after `REPLICATEALL!` because pre-replication
creates scenarios where the dependency matrix has non-zero Float64 entries for
buffer > 0 instances that previously had no dependencies detected.

## 6. Faster Random Benchmark

**File:** `functions/functions_basic.jl`

**What:** Replaced blind random-retry loops in `RANDOMALLOCONCE` and
`RANDOMALLOCMULTI` with direct sampling from valid candidates.

- `RANDOMALLOCONCE`: Instead of randomly picking warehouses until finding one
  with capacity, uses `findall` to get valid warehouses and picks randomly
  from that list.
- `RANDOMALLOCMULTI`: Instead of randomly picking products until finding one
  not yet in the warehouse, builds an `available` list per warehouse and
  samples from it.

**Why:** The original blind-retry approach degenerates when most warehouses are
full (many retries per SKU) or when most products are already placed in a
warehouse (many retries per slot). With 10 warehouses and buffer = 0.2, the
random benchmark took up to 124 seconds. The fix reduces this to seconds.

---

## Benchmark Results (2000 SKUs, all dependency types)

Configuration: `sig = 1e-5`, `min_effect = 0.01`, freq-order REPLICATEALL!,
nor-order Phase 1, WHWEIGHT recomputed after replication.
(Note: the final implementation uses the nor-order ranking for REPLICATEALL!
as well, which performs equal or better per the A/B test described in item 1.)

| Dependency | CHI wins | EMCI wins | Avg diff% |
|------------|-----------|-----------|-----------|
| HD-SF      | 60/60     | 0/60      | -61.3%    |
| HD-VF      | 57/59     | 2/59      | -43.5%    |
| ID-SF      | 38/60     | 22/60     | -0.2%     |
| ID-VF      | 2/59      | 57/59     | +4.1%     |
| MD-SF      | 60/60     | 0/60      | -29.5%    |
| MD-VF      | 47/60     | 13/60     | -10.4%    |
| **Total**  | **264/358** | **94/358** | **73.7%** |

CHI dominates on 5 of 6 dependency types. The remaining weakness on ID-VF
(+4.1%) is a structural floor: with truly independent products and Zipf
frequency, EMCI's pure MCI-based partition is provably optimal (Lin et al. 2025),
and any dependency-aware metric introduces slight suboptimality on this specific
scenario.

## 7. General Local Search (LS) for All Algorithms

**File:** `main_benchmark.jl` (new `RUN_LS!` helper function + calls after each heuristic block)
**Settings:** `main_benchmark_settings.jl` (`apply_ls = true`)

**What:** The local search (`LOCALSEARCHCHI!`) — previously exclusive to CHI — is
now a general post-processing step applied to ALL heuristic algorithms (CHI,
EMCI, GO, GP, GS, BS, KL). When `apply_ls = true`, each algorithm's output W
is copied, improved via pairwise SKU swaps, and recorded as a separate `_LS`
row (e.g., `CHI_1.0e-5_LS`, `EMCI_LS`).

**CHI removed:** The separate CHI algorithm block was removed entirely. CHI was
just CHI + internal local search — now replaced by CHI + general LS, which
produces equivalent results.

**Coappearance matrix:** `Q_ls = COAPPEARENCE(trans_train, sku_weight)` is
computed once per scenario (outside algorithm blocks) and shared by all LS runs.

**Fair timing:** Each `_LS` row's `duration` includes: Q_ls computation time +
original algorithm time + LS execution time. This makes rows self-contained
for comparison.

**LS guard:** Only runs when all SKU weights are uniform
(`all(y->y==sku_weight[1], sku_weight)`), since pairwise swaps assume equal
weight to preserve capacity constraints.

**Key finding:** LS is a powerful equalizer. On HD-SF, EMCI_LS achieves -69%
improvement over base EMCI, nearly closing the gap with CHI. But CHI_LS
still wins by ~5-9% because CHI's initial allocation provides a better
starting point for LS to refine.

## 8. Weight-Aware EMCI (Unified EMCIALLOC)

**File:** `heuristics/heuristic_mci.jl`

**What:** EMCIALLOC now handles both uniform and non-uniform SKU weights in a
single function. It detects weight uniformity internally and branches:

- **Uniform weights:** Uses the original closed-form replication formula from
  Algorithm 5 (Appendix EC.5, Lin et al. 2025). Behavior is identical to the
  previous implementation.
- **Non-uniform weights:** Uses a greedy weight-aware approach:
  1. Ranks SKUs by `MCI / sku_weight[i]` (benefit per unit capacity consumed)
  2. Computes weighted overlap budget: `sum(capacity) - sum(sku_weight)`
  3. Greedily replicates top-ranked SKUs to as many warehouses as the budget
     allows, consuming `sku_weight[i]` per warehouse placement
  4. Assigns remaining SKUs to the warehouse with most remaining capacity

**Why:** The closed-form formula assumes unit weight per SKU and cannot handle
heterogeneous weights. The greedy approach generalizes naturally — when weights
are uniform, it produces functionally equivalent results to the closed form.

**Benefit/weight ranking:** Normalizing MCI by weight ensures heavy SKUs must
provide proportionally more order coverage to justify their warehouse space,
matching CHISQUAREHEUR's approach in `FILLUP!`.

## 9. Variable-Weight Benchmark Mode

**Files:** `main_benchmark.jl`, `main_benchmark_settings.jl`, `run_benchmarks.sh`

**What:** Three SKU weight modes, controlled by `weight_mode` setting:
- `:uniform` — all weights = 1 (default, classic benchmark)
- `:frequency` — weight = number of orders containing that SKU in training data
- `:random` — weight sampled uniformly from `weight_range` (default 1:10)

**Capacity scaling:** For non-uniform weights, capacity is scaled by average
weight to preserve the same fill ratio:
```julia
effective_capacity = round.(Int64, capacity .* avg_weight)
```
A correction step ensures `sum(effective_capacity) >= sum(sku_weight)`.

**Zero-buffer skip:** Constellations with buffer = 0.0 are skipped for
non-uniform weights, since the greedy allocation needs slack to handle weight
variance.

**Algorithm guards:** Algorithms that require uniform weights (KL, IIH, IIHS,
OPT) are automatically skipped when weights are non-uniform. Weight-aware
algorithms (CHI, EMCI, GO, GP, GS, BS, MQKP) run normally.

**Output:** Results CSV includes a `weight_mode` column. Filenames include a
suffix for non-uniform modes (e.g., `2000_benchmark_HD-SF_frequency.csv`).

**Parallel execution:** `run_benchmarks.sh` accepts an optional weight mode
argument. Different weight modes run in separate tmux sessions:
```bash
./run_benchmarks.sh              # uniform (session: benchmarks_2000_uniform)
./run_benchmarks.sh frequency    # (session: benchmarks_2000_frequency)
./run_benchmarks.sh random       # (session: benchmarks_2000_random)
```

## 10. K-Links Stagnation Counter Bug Fix

**Files:** `heuristics/heuristic_klinks.jl`, `functions/functions_k-links.jl`

**What:** Fixed a bug where `STRATEGY1!` and `STRATEGY2!` set `stop = 0` on a
successful move, but since `stop` was passed by value (Julia `Int64`), the reset
never propagated back to the caller's `while stop <= stagnant` loop. The
stagnation counter never actually reset on improvement.

**Fix:** Both strategy functions now **return** the updated stop value (`return 0`
on success, `return stop` otherwise), and the caller captures it:
`stop = STRATEGY1!(X, m, L, capacity, stop)`.

**Impact:** The algorithm now correctly resets the stagnation counter when a
beneficial move is found, matching the intended behavior from Zhang et al. (2021).
Previously, the algorithm always terminated after exactly `stagnant` iterations
regardless of whether improvements were still being found.

## 11. Code Cleanup — Removed Dead Code

**Files:** `functions/functions_chisquare.jl`, `functions/functions_basic.jl`,
`heuristics/heuristic_chisquare.jl`

**Removed unused functions:**
- `ALLOCATE_SMALLCOAPP!` — never called (call was commented out), also contained
  a bug where `allocated[x]` indexed by value instead of position
- `ALLOCATENOCOAPP!` — never called (call was commented out)
- `HYOPTHESISTEST!` — dense-matrix hypothesis test, replaced by
  `HYOPTHESISTEST_SPARSE` which handles both sparse and dense Q
- `DEPTEST!` — only called from `HYOPTHESISTEST!`
- `chi_values!` — only called from `DEPTEST!`
- `fish_values!` — never called (alternative chi test with `independent > 10`
  guard)
- `WHWEIGHT_REVISEDPAPER` — alternative warehouse weight function. A/B tested
  against `WHWEIGHT`: worse by 0.4-1.7% on HD scenarios, noise-level on ID
- `PARCELSSEND_WEIGHT` — never called (2-warehouse-only debug function)

**Removed dense-matrix method dispatches** (Q is always `SparseMatrixCSC` since
`COAPPEARENCE` returns `trans' * trans` on sparse input):
- Dense `REPLICATEALL!` (dep::Matrix, Q::Matrix)
- Dense `ALLOCATEONE!` (dep::Matrix, Q::Matrix)
- Dense `REMOVEALLOC!` (Q::Matrix) — also simplified the call site in
  `heuristic_chisquare.jl` by removing the `if Q isa SparseMatrixCSC` branch
- Dense `CALCVAL` (T::Matrix)
- Dense `REFRESHSTATE!` (Q::Matrix, used `@avxt` from LoopVectorization)

**Other cleanups:**
- Removed unused `freq_order` variable in `heuristic_chisquare.jl`
- Removed trailing orphan comment at end of `functions_chisquare.jl`
- Removed trailing whitespace in `heuristic_greedyseeds.jl`
- Hoisted `nonuniform` weight check out of inner loops in `BESTSELLINGSTART!`
  and `BESTSELLINGTOP!` (was recomputed every iteration)

## 12. Rename CHIM → CHI

**Files:** `main_benchmark.jl`, `main_benchmark_settings.jl`,
`heuristics/heuristic_chisquare.jl`, `helpers/test_benchmark.jl`, `README.md`

**What:** Renamed the chi-square heuristic identifier from `CHIM` to `CHI`
everywhere. The old naming had `CHIM` (without LS) and `CHI` (with built-in LS)
as separate entries. Since local search is now a general post-processing step
applied to all algorithms, the distinction is gone and `CHI` is the natural name.

## 13. Speed Optimizations

**Files:** `functions/functions_chisquare.jl`, `functions/functions_catalan.jl`

**`CURRENTSTATE!`** — replaced per-element `CALCVAL` loop (O(I × K × nnz_per_col))
with a single sparse matrix multiply `state .+= Q * X`. Used in local search
initialization.

**`FILLUP!` state initialization** — same optimization: replaced `CALCVAL` loop
with `Matrix{Float64}(Q * X)`. Also replaced per-iteration `sum(@view(X[i,:]))`
scan with a pre-computed `unallocated` vector. Removed dead dense-Q branch.

**`ADDDEPENDENT!`** — eliminated redundant `argmax`/`findmax` calls (was computing
max of `pot_dep` three times per iteration, now once).

**`REPLICATEALL!` feasibility check** — pre-build sorted `remaining_w` vector
once and maintain incrementally via `searchsortedfirst`/`deleteat!`/`insert!`
instead of allocating a new filtered vector every iteration.

**`BESTSELLINGSTART!` / `BESTSELLINGTOP!`** — same incremental `remaining_w`
optimization as `REPLICATEALL!`.

**`CALCVAL` removed** — no longer called anywhere after `CURRENTSTATE!` and
`FILLUP!` were switched to sparse matrix multiply.

## 14. RAM Optimizations

**Files:** `functions/functions_basic.jl`, `functions/functions_chisquare.jl`

**Coappearance matrix Q: Int64 → Int32** — `COAPPEARENCE` now converts
`trans' * trans` to `SparseMatrixCSC{Int32, Int32}`. Coappearance counts and SKU
indices easily fit in Int32 (max 2.1 billion). Halves Q's memory footprint.

**Dependency matrix dep: Float64 → Float32** — `HYOPTHESISTEST_SPARSE` now builds
the sparse dependency matrix with `Float32` values and `Int32` indices. Dependency
values (`qij - independent`) don't need 15 digits of precision for a heuristic.

**`REMOVEALLOC` rewritten to avoid copying Q** — previously created
`dep_new = convert(SparseMatrixCSC{Float64,Int64}, copy(Q))`, duplicating the
entire Q matrix. Now allocates only a new `Float32` values array and reuses Q's
sparsity structure (`colptr`, `rowval`). Uses merge-scan over Q and dep nonzeros
instead of random sparse lookups.
