# raw_cycle_time_profile

Profile **DIA cycle time across the gradient**, straight from Thermo `.raw` files —
no mzML conversion, no vendor software, scan headers only.

## Why this exists

The nominal duty cycle you plan for is not the one you get, and — more importantly —
**it is often not constant within a run.** A single median cycle time hides that.

When few peptides elute (the start and end of a gradient) the AGC target never fills,
injection time runs to its ceiling, and the cycle stretches. Mid-gradient the ion flux
fills AGC almost immediately, injection time collapses, and the cycle shortens. On a
Stellar DSD sweep this produced a **250% swing inside single runs**:

| retention time | median cycle |
|---|---|
| 0–2 min | 5.21 s (at the MaxIT ceiling) |
| 8–10 min | **1.49 s** |
| 22–24 min | 5.21 s (at the ceiling again) |

A run reported at ~2 s median was actually at 5.2 s during its first and last few
minutes, so anything eluting there got ~4 points across a peak regardless of the summary
number. That is only visible if you keep **every cycle with the retention time it
happened at**, which is what `--per-cycle-out` emits.

## Relationship to `mzml_ion_timing`

Same family of metrics, different input and different granularity. They complement
rather than replace each other:

| | `mzml_ion_timing` | `raw_cycle_time_profile` |
|---|---|---|
| input | mzML | Thermo `.raw` (no conversion step) |
| output | one summary row per file | summary **plus every individual cycle + its RT** |
| cycle detection | precursor-m/z resets | MS1 survey scans, falling back to m/z resets |
| also reports | TIC, ion counts (`injTime × TIC`) | MS2 scan time, injection-time fraction |
| vendor | any | Thermo only |

Use `mzml_ion_timing` if you already have mzML, need ion counts, or are on non-Thermo
data. Use this one to skip conversion, or when you need cycle time resolved against
retention time.

## Dependencies

- `pythonnet`, `pandas`, `numpy` (`pyarrow` only if writing `.parquet`)
- A local **ProteoWizard** install, for its bundled
  `ThermoFisher.CommonCore.Data.dll` and `ThermoFisher.CommonCore.RawFileReader.dll`.
  Auto-detected from the usual install locations; override with `--pwiz <dir>` or
  `$PWIZ_DIR`.

On Windows the DLLs need the .NET Framework runtime — the script sets
`PYTHONNET_RUNTIME=netfx` itself, so no manual setup.

## Usage

```bash
# one summary row per file
python raw_cycle_time_profile.py RAW_DIR --out summary.csv

# plus every individual cycle, for RT-resolved analysis
python raw_cycle_time_profile.py RAW_DIR --out summary.csv \
    --per-cycle-out cycles.parquet

# parse experimental factors out of file names into columns
python raw_cycle_time_profile.py RAW_DIR --out summary.csv \
    --regex '(?P<window>\d+to\d+)_(?P<grad>\d+m)_run(?P<run>\d+)'
```

Accepts any mix of files and directories (recursive by default; `--no-recursive` to
stop that). Duplicate basenames are kept once. `--workers` defaults to 8; a 600-file
sweep takes ~11 min at that setting, and the job is I/O-bound so raising it past your
disk's comfort won't help.

## Output

**Per-file summary** — `File`, `NumScans`, `NumMS1`, `NumMS2`, `ScansPerCycle`,
`CycleDetection`, median/mean/SD/p05/p95/min/max cycle time (s), `MedianMS2ScanTimeMs`,
`MedianMS2InjectionTimeMs`, `InjectionTimeFraction`, `RTStartMin`, `RTEndMin`,
`MS2ScanRange`, `NumCycles`, `Error`, plus any `--regex` named groups.

**Per-cycle** (`--per-cycle-out`) — `File`, `CycleIndex`, `RTMin`, `CycleTimeSec`.
One row per cycle, so ~1,500 rows per 24-min file.

`InjectionTimeFraction` is injection time ÷ per-window scan time — a quick read on
whether a method is **injection-time limited** (approaching 1.0) or **scan-rate
limited** (well below). Useful for knowing which knob will actually move the cycle.

## Caveats

- **Injection times are sampled**, not exhaustive: ~200 MS2 scans per file, reported as
  a median. Per-scan trailer access is the expensive call in this reader, and sampling
  keeps a large sweep tractable. Everything else is computed over all scans.
- **Cycle time is measured, not derived.** It includes real overhead, so it will exceed
  `n_windows × MaxIT`. On Stellar data the per-window overhead measured ~1.7 ms.
- The `mzreset` fallback assumes precursor m/z increases monotonically within a cycle.
  That holds for standard staircase DIA schedules but not for randomized window orders —
  check `CycleDetection` in the output if a file has no MS1 scans.
- Errors are captured per file in the `Error` column rather than aborting the run, so a
  single corrupt file doesn't lose a long sweep. Check that column before trusting counts.

## Provenance

Written for the Stellar DIA DSD Round 2 optimization (July 2026), generalized out of a
one-off analysis of 601 raw files / 582k cycles.
