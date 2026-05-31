# PhylogicNDT Timing Modules: Line-by-Line Walkthrough

This document is a **separate deep-dive explanation** for the classes and methods previously discussed, focusing on:

1. `SinglePatientTiming` execution path and internals
2. `TimingEngine.py` classes and methods used for timing
3. `LeagueModel` execution path and timing aggregation behavior

All file paths below are absolute.

---

## 1) Entry point wiring (`PhylogicNDT.py`)

File: `/tmp/workspace/shwong-tw/PhylogicNDT/PhylogicNDT.py`

### Lines 405–423: timing CLI registration
- **405–406**: Adds `Timing` subparser.
- **407–412**: Defines `-min_supporting_muts` argument used by timing.
- **413**: `timing.set_defaults(func=SinglePatientTiming.SinglePatientTiming.run_tool)` binds command to the timing module.
- **415–423**: Adds alias command `SinglePatientTiming` with same `min_supporting_muts`, also mapped to `run_tool`.

### Lines 425–505: league CLI registration
- **425–426**: Adds `LeagueModel` subparser.
- **427–504**: Defines cohort/comparison/permutation options.
- **505**: `leaguemodel.set_defaults(func=LeagueModel.LeagueModel.run_league_model)` binds command to league module.

### Lines 513–515: dispatch
- **513**: Parse args.
- **515**: Execute the function stored in `args.func`.

---

## 2) SinglePatientTiming module

File: `/tmp/workspace/shwong-tw/PhylogicNDT/SinglePatientTiming/SinglePatientTiming.py`

## Method: `run_tool(args)` (lines 9–58)

- **10**: Logs parsed arguments.
- **12–14**: Constructs `Patient` object using patient ID and blacklist/whitelist/driver files.
- **16**: Branch: if `--sif` provided.
- **17–28**: Reads each row of sample info file, parses `sample_id/maf_fn/seg_fn/purity/timepoint`, and calls `patient_data.addSample(...)` with `input_type='post-clustering'` and `seg_input_type='timing_format'`.
- **30**: Else branch uses `-s sample_id:maf:seg:purity:timepoint` entries.
- **34–45**: Splits each entry and calls `addSample(...)` similarly.
- **47**: `patient_data.preprocess_samples()` synchronizes variants across samples.
- **48–49**: Validates `min_supporting_muts >= 1`.
- **50**: Creates `TimingEngine` with patient and threshold.
- **51**: Calls `timing_engine.time_events()` to compute event timing distributions.
- **52–53**: Writes main timing output TSV.
- **54–55**: Loads driver genes.
- **56**: Builds pairwise event-order probabilities via `compare_events`.
- **57**: Writes `.comp.tsv` used by LeagueModel.

## Method: `compare_events(timing_engine, drivers=())` (lines 61–113)

- **65**: Initialize event list.
- **66–68**: Include WGD if present.
- **68**: Include all CN events across regions.
- **69–71**: Include driver coding mutations (protein-changing heuristic).
- **73**: Iterate all event pairs.
- **74–75**: Skip pairs missing timing distributions.
- **76**: Compute difference between two pi distributions.
- **77–102**: Integrate signed difference from both ends to estimate probability that event1 is before event2, event2 before event1, or unknown/overlap.
- **112**: Save tuple `(p_before, p_after, p_unknown)` keyed by `(event1,event2)`.
- **113**: Return pairwise comparison dictionary.

---

## 3) TimingEngine classes and methods

File: `/tmp/workspace/shwong-tw/PhylogicNDT/SinglePatientTiming/TimingEngine.py`

## Globals (lines 11–16)
- Define chromosomes (`1..22`,`X`), arm labels (`p`,`q`), default CN whitelist (`(1,2),(0,2),(2,2)`), chromosome sizes and centromere map.

## Class: `TimingEngine` (starts line 19)

### `__init__` (23–43)
- Saves patient/config.
- Builds `TimingSample` for each sample.
- Initializes stores (`concordant_cn_states`, `timeable_muts`, `all_cn_events`, etc.).
- Runs pipeline in constructor:
  1. `get_concordant_cn_states()`
  2. `call_wgd()`
  3. `get_mutations()`
  4. `get_arm_level_cn_events()`

### `get_concordant_cn_states` (44–100)
- For each arm region in each sample:
  - Coerces CNs near integers to discrete states.
  - Keeps only states present in **all samples**.
  - Intersects supporting mutation IDs across samples.
  - Rebuilds each shared mutation as multisample `TimingMut` (arrays for alt/ref/CN/CCF across samples).
  - Creates multisample `TimingCNState` per concordant arm.

### `call_wgd` (101–128)
- For each sample:
  - Rebuild concordant mutation interval trees.
  - Recompute arm states on concordant set.
  - Call sample-level WGD and concordant arm events.
  - If sample has concordant WGD, collect sample WGD pi distribution.
- Requires all samples to have concordant WGD and minimum concordance threshold.
- Builds patient-level `TimingWGD` from supporting regions with major copy >=2 in all samples.
- Re-calls CN-state events in WGD-relative mode (`baseline=2`).

### `get_mutations` (129–150)
- Starts from sample 0 lookup order.
- Reuses already timeable concordant mutations when possible.
- Otherwise rebuilds multisample `TimingMut` from per-sample arrays.

### `get_arm_level_cn_events` (152–176)
- Initializes `truncal_cn_events` and `all_cn_events` keyed by `gain_/loss_ + chrom + arm`.
- Uses `_get_cluster_ccfs()` to infer clonal concordance of each CN event against cluster CCF profiles.
- Marks each CN event `is_clonal`.
- Populates all events and truncal-only subset.

### `_get_cluster_ccfs` (177–187)
- Aggregates log CCF densities by mutation cluster ID.
- Normalizes each sample’s cluster profile.
- Returns dictionary `cluster_id -> [n_samples x 101]`.

### `time_events` (189–228)
- Defines defaults:
  - uniform distribution for unknown clonal timing.
  - delta-at-1 for subclonal events.
- If WGD exists:
  - Time WGD first.
  - Time clonal mutations relative to WGD.
  - CN events: subclonal gets delta-at-1; clonal gains/losses use WGD-based methods.
- If no WGD:
  - Subclonal CN events get delta-at-1.
  - Clonal gains with enough support and whitelisted state use gain model.
  - Supporting clonal mutations are timed relative to that gain.
  - Otherwise CN gets uniform timing.
- Any mutation still unassigned gets uniform (if clonal) or subclonal delta.

---

## Class: `TimingSample` (starts line 230)

### `__init__` (234–259)
- Stores sample metadata and CN profile.
- Initializes mutation interval trees and lookup tables.
- Calls:
  1. `fill_mutation_intervaltree()`
  2. `call_arm_level_cn_states()`
  3. `call_wgd()`
  4. `get_arm_level_cn_events()`

### `fill_mutation_intervaltree` (266–297)
- Iterates high-confidence + low-coverage mutations.
- Retrieves local CN from segment tree at mutation position.
- Creates `TimingMut` per mutation.
- Computes multiplicity distribution when local CN available.
- Stores in sample-wide interval tree/lookup; or in concordant-only tree if `concordant_muts` passed.

### `call_arm_level_cn_states` (298–352)
- For each chrom arm:
  - Slices segment tree by centromere (`p` vs `q`).
  - Computes bp coverage and bp per CN state.
  - Calls arm state if one state dominates >40% observed arm bp.
  - Marks arm missing if <50% of true arm covered.
  - Creates `TimingCNState` with supporting mutations matching that state.

### `call_wgd` (353–386)
- `use_concordant_states=True` path:
  - Requires pre-existing sample WGD.
  - Collects concordant states supporting WGD (0/2 or >=2/2 logic).
  - Creates `concordant_WGD` and calls concordant events in WGD mode.
- Default path:
  - Uses sample arm states to count WGD support.
  - Calls WGD if enough both-allele gains and broad support across evaluable arms.
  - Re-calls state events relative to WGD baseline.

### `get_arm_level_cn_events` (387–400)
- Resets event dictionary.
- Flattens each state’s `cn_events` into map keyed by event name.

---

## Class: `TimingWGD` (starts line 403)

### `__init__` (410–416)
- Stores supporting arm states and parameters.
- Builds combined list of supporting mutations.

### `get_pi_dist` (417–433)
- For each supporting arm gain event with enough mutations:
  - Ensure gain pi distribution is computed.
  - Add log-probabilities across events.
- Normalize to produce WGD pi distribution.

### `get_supporting_muts` (434–438)
- Concatenates supporting mutations from all supporting states.

---

## Class: `TimingCNState` (starts line 440)

### `__init__` (445–454)
- Stores region and allelic CN values.
- Immediately calls `call_events(...)`.

### `p2_domain` property (458–465)
- Defines valid p2 domain for specific CN states where model is implemented.

### `call_events` (466–496)
- Sets baseline CN to 1 (no WGD) or 2 (WGD-relative).
- Rounds/coerces allele CN around baseline.
- If allele1 differs from baseline, emits gain/loss `TimingCNEvent` with `ccf_hat`.
- Same for allele2.
- This can create 0, 1, or 2 events per state.

---

## Class: `TimingCNEvent` (starts line 498)

### `__init__` (499–518)
- Stores event metadata and support.
- Parses copy number string into numeric `cn_a1/cn_a2`, stores `allelic_cn`, `total_cn`.

### `event_name` property (523–527)
- Converts internal type like `Arm_gain` to output key like `gain_8p`.
- Raises `NotImplementedError` for non-arm type prefixes.

### `log_p2_prior` property (530–536)
- Prior over p2 for implemented gain states.

### `get_p2_dist_for_gain` (538–552)
- Combines mutation multiplicity evidence across supporting clonal mutations.
- Applies prior and normalization.
- Corrects p2 using simulation-based transform.

### `get_pi_dist_for_gain` (554–593)
- Rejects non-whitelisted states.
- Converts corrected p2 density to pi-space by change of variables.
- Separate transform paths for 0/2-or-2/2-like vs 1/2-like states.
- Clips/norms final pi distribution.

### `_correct_p2` (594–617)
- Uses simulated detection distortion curve to map posterior p2 to corrected “real” p2 density.

### `_simulate_p2` (618–650)
- Simulates mutations under varying p2, with sample purity and CN-dependent expected AF.
- Applies detectability criteria to estimate observed p2 mapping.

### `get_pi_dist_for_loss` (652–666)
- Uses WGD CDF to set loss timing:
  - allelic CN 0 -> before WGD
  - allelic CN 1 -> after WGD

### `get_pi_dist_for_higher_gain` (667–681)
- Uses WGD CDF for gains above baseline-2:
  - allelic CN 3 -> after WGD
  - allelic CN >=4 -> before WGD

### `log_timing_info` (682–684)
- Appends unique timing messages.

---

## Class: `TimingMut` (starts line 687)

### `__init__` (691–718)
- Stores mutation metadata and per-sample arrays.
- Computes `ccf_hat` and clonal flag:
  - cluster assignment 1 implies clonal;
  - otherwise clonal by CCF threshold.
- Initializes multiplicity likelihood containers.

### `var_str` property (730–733)
- Canonical mutation key `chr:pos:ref:alt`.

### `event_name` property (735–738)
- Prefer `GENE_PROTEINCHANGE`, else `var_str`.

### `get_multiplicity_likelihoods` (740–752)
- For multiplicity values 1..max major CN, computes binomial likelihood of observed alt counts.

### `get_log_mult_1_2_distribution` (754–759)
- Builds log distribution interpolating multiplicity 1 vs 2.

### `get_pi_dist` (761–796)
- Times mutation relative to a matched gain event.
- Requires stable local CN state and multiplicity likelihoods.
- Computes before/after-gain likelihood and combines with gain CDF into mutation pi distribution.

### `get_mult_dist` (797–807)
- Validates consistent local CN.
- Computes multiplicity likelihood dictionary.
- For CN states `(0|1|2, 2)`, caches 1-vs-2 multiplicity log distribution.

---

## 4) Output layer for timing

File: `/tmp/workspace/shwong-tw/PhylogicNDT/output/PhylogicOutput.py`

## `write_timing_tsv(self, timing_engine)` (1129–1189)

- **1135–1139**: Builds header and patient name.
- **1141**: Opens `<patient>.timing.tsv`.
- **1143–1147**: Writes WGD row first if present.
- **1148–1159**: Writes each CN event row, with pi summary and full pi bins.
- **1160–1170**: Writes mutation rows similarly.
- **1171–1189**: If WGD exists, writes separate `<patient>.WGD_supporting_events.timing.tsv`.

## `write_comp_table(indiv_id, comps)` (1191–1197)
- Writes `<indiv>.comp.tsv` rows: sample, event pair, p(event1 before 2), p(event2 before 1), unknown.

---

## 5) LeagueModel module

File: `/tmp/workspace/shwong-tw/PhylogicNDT/LeagueModel/LeagueModel.py`

## `run_league_model(args)` (3–140)
- **10–25**: Chooses input source:
  - aggregated comparison file (`--comparison_fn`) or
  - concatenated list of `.comp.tsv` files (`--comps`).
- **27**: Reads comparisons table to DataFrame.
- **28–33**: Optional final sample list parsing.
- **34–38**: Instantiates `LeagueModelData.League(...)`.
- **39–44**: Writes `<cohort>.all_events.tsv` (full per-sample events).
- **46–48**: Full-run baseline + initialize odds store.
- **51–76**: Writes matrix of pairwise comparison probabilities.
- **77–84**: For each permutation:
  - random sample subset,
  - run league simulation,
  - update odds.
- **85**: Compute final log-odds arrays.
- **88–95**: Plot and save odds/position figures.
- **98–120**: Write prevalence outputs.
- **122–138**: Write log-odds table.
- **140**: Pickle full league object.

---

## 6) LeagueModel data structures and algorithms

File: `/tmp/workspace/shwong-tw/PhylogicNDT/LeagueModel/LeagueModelData.py`

## Class: `Eve_Pair` (15–38)
- Stores an event pair and win-rate counts (`event1`, `event2`, `unknown`).
- `calculate_rates()` converts counts to probabilities.

## Class: `Season` (39–84)
- Holds season points table across events.
- `update_table_from_multinomial_sampling`: awards 2/1/0 points depending on draw outcome.
- `get_sorted_league`: rank by points.
- `get_event_positions`: assign rank positions with tie handling.

## Class: `League` (85 onward)

### Static `arm_CNVs` construction (88–97)
- Predefines arm-level gain/loss labels and includes `WGD`.

### `__init__` (98–137)
- Stores inputs/options.
- Calls pipeline:
  1. `load_df()`
  2. `subset_to_earliest_point_mut()`
  3. optional sample filtering
  4. `update_pairwise_probs()`
  5. `calc_event_occur()` and `get_final_event_list()` (unless provided)
  6. `form_pairs_for_league_model()`
  7. `update_pairs_for_league_model()`
  8. `run_league_model_iter(num_seasons=1000)`

### `load_df` (218–271)
- Reads each row (sample, event1, event2, probabilities).
- Creates per-sample event sets and pairwise probability maps.
- Collapses event names to gene/event-level keys.
- Classifies event type into arm-level/focal/WGD/snv.

### `subset_to_earliest_point_mut` (276–305)
- Keeps all CN/WGD/focal events.
- For multiple point mutations in same gene, keeps earliest-supported hit by comparison summary.

### `update_pairwise_probs` (306–325)
- Rebuilds per-sample pair probabilities using filtered final events.

### `calc_event_occur` (144–162)
- Counts event prevalence over selected samples.

### `get_final_event_list` (326–387)
- Builds final panel in priority order:
  1. common arm gains
  2. common arm losses
  3. remaining arm events
  4. SNV events
  5. focal gains/losses
  6. homdel
  7. WGD
- Applies prevalence thresholds and caps.

### `form_pairs_for_league_model` (389–399)
- Creates all event pairs in final panel.

### `update_pairs_for_league_model` (400–425)
- Aggregates per-sample pairwise probabilities into global pair objects.
- Adds pseudo-count behavior when support is low.
- Finalizes pair win rates.

### `run_league_model_iter` (426–447)
- Repeatedly simulates seasons.
- For each pair, draws multinomial outcome from pair win/draw rates.
- Updates rankings and records event positions.

### `calc_odds` (449–470)
- Converts event position distributions into early/late odds.

### `init_odds` (472–475)
- Initializes odds history arrays for each event.

### `run_permutation` (476–485)
- Resets season state and reruns league simulation for one subset/permutation.

### `update_odds` (486–493)
- Appends current permutation odds to event histories.

### `run_full_run` (494–503)
- Executes baseline run and stores `odds_full_run` and full prevalence.

### `calc_log_odds_full_run` (504–509)
- Takes log10 of early-odds history per event.

### `plot_league_run` (522 onward)
- Produces odds or positional plots and color-codes event types.

---

## 7) Short interpretation summary

- `SinglePatientTiming` computes per-patient timing distributions (`pi_dist`) for WGD, arm-level CN events, and mutations.
- `compare_events` converts distributions to pairwise order probabilities.
- `LeagueModel` reads those pairwise probabilities across samples and uses simulation-based ranking to infer cohort-level relative ordering.
- CN-event handling in LeagueModel is label/probability aggregation; mechanistic timing is done upstream by `TimingEngine`.
