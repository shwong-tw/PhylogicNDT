# Focal Copy Number Region Support – Detailed Modification Guide

This document describes **every location** in the codebase that was modified to add support for focal (sub-arm-level) copy number events driven by an external peaks file (e.g., GISTIC output).

---

## Table of Contents

1. [Overview](#overview)
2. [File: `PhylogicNDT.py`](#file-phylogicndtpy)
3. [File: `Cluster/Cluster.py`](#file-clusterclusterpy)
4. [File: `SinglePatientTiming/SinglePatientTiming.py`](#file-singlepatienttimingsinglepatienttimingpy)
5. [File: `SinglePatientTiming/TimingEngine.py`](#file-singlepatienttimingtimingengine-py)
6. [File: `data/Patient.py`](#file-datapatientpy)
7. [File: `data/SomaticEvents.py`](#file-datasomaticevents-py)

---

## Overview

The goal of these changes is to allow users to supply a file describing **focal CN regions** (gains, losses, homozygous deletions) via the `--cn_peaks` CLI argument. The system then:

1. Parses the file (tab- or comma-delimited).
2. Computes weighted-average CN across the region for each sample.
3. Creates `CopyNumberEvent` objects tagged as `Focal_gain`, `Focal_loss`, or `Focal_homdel`.
4. Integrates these events into the timing engine so they receive pi-distributions and are timed alongside arm-level events.

---

## File: `PhylogicNDT.py`

### Location: Lines ~412–417 (inside the `timing` subparser block)

```python
timing.add_argument('--cn_peaks',
                    type=str,
                    action='store',
                    dest='gistic_fn',
                    default=None,
                    help='Tab/comma-delimited focal region file with columns event_class, chromosome, start, end (optional region_label)')
```

**Purpose:** Adds the `--cn_peaks` command-line argument to the **`Timing`** subcommand so users can pass a focal regions file.

---

### Location: Lines ~428–433 (inside the `single_patient_timing` subparser block)

```python
single_patient_timing.add_argument('--cn_peaks',
                    type=str,
                    action='store',
                    dest='gistic_fn',
                    default=None,
                    help='Tab/comma-delimited focal region file with columns event_class, chromosome, start, end (optional region_label)')
```

**Purpose:** Same argument added to the **`SinglePatientTiming`** subcommand (alias).

---

## File: `Cluster/Cluster.py`

### Location: Line ~79 (after `patient_data.get_arm_level_cn_events()`)

```python
if args.gistic_fn:
    patient_data.get_focal_level_cn_events(args.gistic_fn)
```

**Purpose:** When a focal regions file is supplied during clustering, focal CN events are extracted from the patient's segment data before pre-processing. This ensures that focal events appear in the clustering results alongside arm-level events.

---

## File: `SinglePatientTiming/SinglePatientTiming.py`

### Location: Line ~50 (TimingEngine instantiation)

**Before:**
```python
timing_engine = TimingEngine.TimingEngine(patient_data, min_supporting_muts=args.min_supporting_muts)
```

**After:**
```python
timing_engine = TimingEngine.TimingEngine(patient_data, min_supporting_muts=args.min_supporting_muts,
                                          focal_cn_fn=args.gistic_fn)
```

**Purpose:** Passes the focal CN filename through to the `TimingEngine` constructor so that focal events can be loaded and timed.

---

## File: `SinglePatientTiming/TimingEngine.py`

### Location: Line 3 – New import

```python
import csv
```

**Purpose:** Needed for parsing the focal CN peaks file.

---

### Location: Line ~24 – Constructor signature change

**Before:**
```python
def __init__(self, patient, cn_state_whitelist=_cn_state_whitelist, chromosomes=_chromosomes,
             arms=_arms, min_supporting_muts=3):
```

**After:**
```python
def __init__(self, patient, cn_state_whitelist=_cn_state_whitelist, chromosomes=_chromosomes,
             arms=_arms, min_supporting_muts=3, focal_cn_fn=None):
```

**Purpose:** Accepts the optional focal CN filename parameter.

---

### Location: Lines ~42–43 – Conditional call in `__init__`

```python
if focal_cn_fn:
    self.get_focal_level_cn_events(focal_cn_fn)
```

**Purpose:** If a focal CN file was provided, load and register those events immediately after arm-level events are loaded.

---

### Location: Lines ~179–267 – New static methods and `get_focal_level_cn_events`

#### `_normalize_chromosome(chromosome)` (static, line ~179)

Strips `chr` prefix, converts `23`→`X`, `24`→`Y`. Ensures consistent chromosome naming.

#### `_normalize_focal_event_class(event_class)` (static, line ~189)

Maps user-supplied event classes (e.g., `amp`, `del`, `deepdel`) to canonical internal names (`gain`, `loss`, `homdel`).

#### `_load_focal_regions(self, cn_peaks_fn)` (line ~207)

- Opens the peaks file.
- Auto-detects tab vs comma delimiter.
- Resolves column names through flexible alias matching (e.g., `chr` → `chromosome`).
- Validates required columns and parses rows into a list of region dictionaries.

#### `get_focal_level_cn_events(self, cn_peaks_fn)` (line ~266)

- Iterates over loaded regions.
- For each region, queries each sample's `CnProfile` interval tree for overlapping segments.
- Computes bp-weighted average CN for allele 1 and allele 2.
- Determines whether the CN event is consistent with a gain/loss/homdel across all samples.
- Determines clonality by comparing CCF to cluster CCF distributions.
- Creates a `TimingCNEvent` and inserts it into `self.all_cn_events` and (if clonal) `self.truncal_cn_events`.

---

### Location: Lines ~377–379, ~389–391 – Pi-distribution assignment for focal events

Inside `time_events()` (two code paths: with WGD and without WGD):

```python
if cn_event.Type.startswith('Focal_'):
    cn_event.pi_dist = uniform_dist if cn_event.is_clonal else subclonal_dist
    continue
```

**Purpose:** Focal events don't have supporting mutations to estimate timing from (they span arbitrary small regions), so they receive a **uniform** prior if clonal, or a **subclonal** prior otherwise. The `continue` prevents them from falling into the arm-level gain/loss timing logic which relies on supporting mutation counts.

---

### Location: Lines ~709–710 – `TimingCNEvent.event_name` property

```python
if self.Type.startswith('Focal_'):
    return self.Type[6:] + '_' + self.arm
```

**Purpose:** Generates human-readable event names for focal events (e.g., `gain_MYC:8:12345-67890`). The `self.arm` field stores the region label for focal events. The `[6:]` strips the `Focal_` prefix to yield `gain`, `loss`, or `homdel`.

---

## File: `data/Patient.py`

### Location: Line ~12 – New import

```python
import csv
```

**Purpose:** Parsing the focal CN peaks file within the Patient class.

---

### Location: Lines ~369–404 – Static helper methods

#### `_normalize_chromosome(chromosome)` (line ~369)

Identical logic as in `TimingEngine`; strips `chr`, converts numeric sex chromosomes.

#### `_normalize_focal_event_class(event_class)` (line ~380)

Identical logic as in `TimingEngine`; maps aliases to canonical class names.

---

### Location: Lines ~406–461 – `_load_focal_regions(self, cn_peaks_fn)`

Same parsing logic as in `TimingEngine`:
- Auto-detects delimiter.
- Resolves column aliases.
- Validates rows and returns a list of region dicts.

Used by the Patient-level focal loading path (invoked from `Cluster.py`).

---

### Location: Lines ~463–558 – `get_focal_level_cn_events(self, cn_peaks_fn)`

Patient-level method (used during clustering):
- Iterates regions.
- For each sample, queries `CnProfile[chrom][start:end]` for overlapping segments.
- Computes bp-weighted average CN and CCF values per allele.
- Validates that ≥50% of the region is covered (`overlap_bp < region_len * 0.5` → skip).
- Selects the appropriate allele based on event class (max for gain, min for loss/homdel).
- Calls `self._add_cn_event_to_samples(...)` with a `region_label` keyword argument.

---

### Location: Line ~816 – `_add_cn_event_to_samples` signature change

**Before:**
```python
def _add_cn_event_to_samples(self, chrom, start, end, arm, cns, cn_category, ccf_hat, ccf_high, ccf_low):
```

**After:**
```python
def _add_cn_event_to_samples(self, chrom, start, end, arm, cns, cn_category, ccf_hat, ccf_high, ccf_low,
                             region_label=None):
```

**Purpose:** Passes an optional `region_label` through to the `CopyNumberEvent` constructor.

---

### Location: Line ~828 – Passing `region_label` to `CopyNumberEvent`

**Before:**
```python
cn = CopyNumberEvent(..., arm=arm)
```

**After:**
```python
cn = CopyNumberEvent(..., arm=arm, region_label=region_label)
```

**Purpose:** Propagates the region label to the event object so it can be used in naming and output.

---

## File: `data/SomaticEvents.py`

### Location: Line ~243 – Constructor signature change

**Before:**
```python
def __init__(self, chrN, cn_category, ..., dupe=False):
```

**After:**
```python
def __init__(self, chrN, cn_category, ..., dupe=False, region_label=None):
```

**Purpose:** Accepts an optional `region_label` parameter.

---

### Location: Line ~255 – Store `region_label` attribute

```python
self.region_label = region_label
```

**Purpose:** Persists the label on the event object for downstream use in `event_name`.

---

### Location: Lines ~288–290 – `event_name` logic for focal events

**Before:**
```python
elif cn_category.startswith('Focal'):
    gl = cn_category.split('_')[1]
    self.event_name = gl + '_' + str(chrN) + start.band + '-' + end.band[1:] if start != end else gl + str(chrN) + start.band
```

**After:**
```python
elif cn_category.startswith('Focal'):
    gl = cn_category.split('_')[1]
    if self.region_label is None:
        self.region_label = '{}:{}-{}'.format(self.chrN, self.start, self.end)
    self.event_name = gl + '_' + str(self.region_label)
```

**Purpose:** The old code assumed `start`/`end` had `.band` attributes (cytoband objects). The new code uses `region_label` (either user-supplied or auto-generated from coordinates) to form the event name, avoiding errors with integer coordinates.

---

## Summary Table

| File | Line(s) | Change Type | Description |
|------|---------|-------------|-------------|
| `PhylogicNDT.py` | ~412–417, ~428–433 | CLI argument | Adds `--cn_peaks` to Timing and SinglePatientTiming subcommands |
| `Cluster/Cluster.py` | ~79–80 | Conditional call | Calls `get_focal_level_cn_events` during clustering |
| `SinglePatientTiming/SinglePatientTiming.py` | ~50 | Param passing | Forwards `focal_cn_fn` to `TimingEngine` |
| `SinglePatientTiming/TimingEngine.py` | 3 | Import | Adds `csv` import |
| `SinglePatientTiming/TimingEngine.py` | ~24 | Signature | Adds `focal_cn_fn` param to constructor |
| `SinglePatientTiming/TimingEngine.py` | ~42–43 | Init logic | Conditionally loads focal events |
| `SinglePatientTiming/TimingEngine.py` | ~179–267 | New methods | `_normalize_chromosome`, `_normalize_focal_event_class`, `_load_focal_regions`, `get_focal_level_cn_events` |
| `SinglePatientTiming/TimingEngine.py` | ~377–379, ~389–391 | Pi-dist | Assigns uniform/subclonal dist to focal events |
| `SinglePatientTiming/TimingEngine.py` | ~709–710 | Property | `event_name` for focal `TimingCNEvent` |
| `data/Patient.py` | ~12 | Import | Adds `csv` import |
| `data/Patient.py` | ~369–404 | New methods | `_normalize_chromosome`, `_normalize_focal_event_class` |
| `data/Patient.py` | ~406–461 | New method | `_load_focal_regions` |
| `data/Patient.py` | ~463–558 | New method | `get_focal_level_cn_events` |
| `data/Patient.py` | ~816 | Signature | Adds `region_label` to `_add_cn_event_to_samples` |
| `data/Patient.py` | ~828 | Param passing | Passes `region_label` to `CopyNumberEvent` |
| `data/SomaticEvents.py` | ~243 | Signature | Adds `region_label` to `CopyNumberEvent.__init__` |
| `data/SomaticEvents.py` | ~255 | Attribute | Stores `self.region_label` |
| `data/SomaticEvents.py` | ~288–290 | Event naming | Uses `region_label` for focal event names |
