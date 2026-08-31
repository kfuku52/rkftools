# Performance

Measurements below compare commit `81d94e5` (0.1.10) with the 0.1.11 changes on
the same macOS arm64 / R 4.4.3 environment, using ape 5.8-1 and phytools 2.5-2.
Each case uses deterministic inputs, one untimed warm-up, and three measured
runs. Values are medians. These are local measurements, not universal timing
promises; shared-machine load and garbage collection affect absolute times.

| Case | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| Table conversion, 800-tip comb | 15.707 s | 0.358 s | 44x |
| Table conversion, 800-tip random | 0.847 s | 0.233 s | 3.6x |
| Root mapping, 800-tip random | 1.812 s | 0.017 s | 107x |
| Root mapping, 800-tip comb | 1.005 s | 0.011 s | 91x |
| MAD, 24-tip random | 0.135 s | 0.064 s | 2.1x |

All 19 benchmark outputs matched the baseline within `1e-8`, including complete
branch tables and MAD custom results. This equivalence comparison uses inputs
where the original implementation was correct; dedicated regression tests
cover intentional corrections to padding, missing observations, and duplicate
MAD tips.

The complete benchmark process, including package loading, correctness checks,
warm-up, garbage collection, and repetitions, took 112.43 s before and 58.37 s
after. `/usr/bin/time -l` measured maximum process RSS of 503.6 MiB before and
376.3 MiB after. Those are whole-process peaks, not per-function allocations.
Per-case CSV files also report the sum of R's two maximum-used heap columns
from `gc()`. That heap estimate excludes native allocations and is distinct
from process RSS.

The table reports the final full run, including its shared-host timing noise.
Earlier complete runs measured the 800-tip comb table conversion at 0.107 and
0.140 s, and random-tree root mapping at 0.008 and 0.013 s. All runs preserved
the benchmark outputs. The final run's 100-tip star scoring case triggered the
2x gate (0.053 to 0.180 s). Nine alternating before/after measurements did not
reproduce that regression: medians were 0.207/0.157 s wall time and
0.068/0.072 s user CPU. The corresponding 200-tip control measured
0.462/0.550 s wall time and 0.253/0.283 s user CPU. These controls support
host-load sensitivity rather than a threefold code slowdown; they do not claim
that every operation became faster. The gate remains unchanged.

Root split matching now counts descendant tips and target-side membership in
one traversal. Table conversion reuses descendants and preallocates columns;
its historical clade-bit ordering still requires quadratic work for very large
trees. The change removes repeated whole-edge scans, not every possible scaling
limit. MAD no longer forces a full garbage collection on every serial call.

## Reproduce or compare a change

`make benchmark` always loads the current checkout with pkgload, verifies the
namespace source path, and logs the package version and commit. It cannot
silently benchmark an older globally installed package.

```sh
make setup-minimal
make benchmark

# Compare against the previous implementation without switching branches.
benchmark_before=$(mktemp -d)
git archive 81d94e5 | tar -x -C "$benchmark_before"
RKFTOOLS_BENCHMARK_COMMIT=81d94e5 Rscript tools/benchmark.R \
  "--source=$benchmark_before" --output=benchmark/before.csv
Rscript tools/benchmark.R --output=benchmark/after.csv \
  --baseline=benchmark/before.csv
```

Each CSV has a companion `.samples.csv`, `.outputs.rds`, and `.metadata.rds`.
The output file stores complete typed results for equivalence checks; keep it
with the baseline CSV. A comparison also writes `.comparison.csv`.

Use `--sizes=200,400,800`, `--repetitions=5`, or an output path to customize a
run. `--max-ratio=2.0` makes a baseline comparison fail if a median more than
doubles and the new median exceeds 50 ms, avoiding ratios dominated by timer
resolution. Compare matching environments and repeat a flagged case before
attributing it to a code change. Different backend versions can also change
correct results and require explicit review of a new baseline.

The weekly/manual benchmark workflow retrieves the most recent successful
run's retained artifact as its baseline. The first run (or an expired artifact)
establishes a baseline. It checks output equivalence and uses the conservative
2x/50 ms regression gate, with a configurable manual ratio. CSV, typed outputs,
metadata, and process RSS logs are retained for fourteen days. Routine PR tests
continue to enforce correctness independently of benchmark timing noise.
