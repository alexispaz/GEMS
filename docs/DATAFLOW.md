# GEMS Dataflow

Related documents: [architecture](ARCHITECTURE.md), [memory model](MEMORY_MODEL.md), [input language](INPUT_LANGUAGE.md), and [CLI reference](CLI_REFERENCE.md).

## Startup

`Main.F90` initializes:

- random number generator and elements;
- output precision;
- parser variables (`jobname`, `i`, `ci`, `boxx`, `boxy`, `boxz`, `time`, `step`, `dt`, `pi`);
- group registry `gindex`;
- calculation registry `cindex`;
- `sys`, `ghost`, `gsel`, `gmeta`, and saved groups `gr(1:mgr)`;
- interaction registry `ngindex`, output registry `of_vop`, and integration registry `its`;
- parser options, then the input read/execute loop.

## Input To State

```text
read_line
    -> directives / variables / math expansion / blocks
    -> first token
    -> execute_command
    -> group/atom/interactions/integrators/output mutations
```

Creation commands (`+`, `>+`, `^+`) allocate atoms through selection/create helpers and attach them to `sys` and/or `gsel`.

Selection commands (`>`, `^`, `^>`, `^~`) manipulate `gsel` through group attach/detach.

Saved groups are created by `group i add`.

## Interaction Setup

`interact` command flow:

```text
interact [:label] ref_group [< b_group] kind subkind params...
    -> interact_new
    -> allocate concrete ngroup extension
    -> init ngroup
    -> attach gr(ref) to g%ref
    -> attach gr(b) to g%b
    -> attach gr(ref/b) to ngroup index
    -> link label
    -> concrete cli
    -> test_update
```

Examples:

```gms
interact 1 pair lj 1.0 1.0 2.5
interact 2 < 3 pair lj 0.1 3.4 200
interact :gr 1 graph subgraphs 3.0
```

Unlike the older architecture, this version does not show a global `nomoreigs` flag. Interaction setup immediately calls `test_update`, so it can build/touch neighbor structures during setup.

## Dynamics Loop

`Programs.f90::dinamic`:

```text
interact(.false.)

for each step:
    dm_steps += 1
    integration_stepa()
    interact(b_out)
    integration_stepb()
    optional write_out common phase
    optional checkpoint
    optional scheduled block
    optional timing
```

Integration step A/B dispatches procedure pointers stored in `its`. Both steps call `gindex_all_changed`.

## Interaction Runtime

`Interaction.f90::interact`:

```text
test_update()
zero bias/energy
zero forces and epot on sys and ghost unless dontclean
for each ngroup in ngindex:
    call g%interact()
    accumulate energy
for each force-field bound group:
    call interact()
optional write_out phase 2
CV calc hooks
```

## Neighbor Update

`test_update`:

1. Applies PBC/ghost updates: `pbchalfghost(maxrcut)` when `useghost`, otherwise `do_pbc(sys)`.
2. Tessellates and sorts each ngroup candidate group `b`.
3. If any ngroup is not listed, calls `update`.
4. Otherwise computes displacement maxima and calls `update` when the Verlet skin threshold is exceeded.

`update`:

1. Increments `nupd_vlist`.
2. Saves `pos_old` for all atoms in `sys`.
3. Selects each ngroup neighbor algorithm (`cells` or `verlet`) based on candidate group tessellation.
4. Builds `nn/list`.

## Output And Analysis

`outfile` commands create persistent output objects. `write_out` runs during common, interaction, force-field, and specialized phases.

`out state` is immediate: it calls `interact(.true.)` and writes force output. This can build/update neighbor structures as a side effect.

CVs are configured with `cv` commands and are evaluated after interactions.

## Lifecycle Boundaries

Important practical boundaries:

- Atom creation/selection: linked group state changes.
- Interaction declaration: ngroups are allocated and indexed.
- First `test_update`: ghost, cell, and neighbor structures can be built.
- Each integration step: positions/velocities/box may change and group caches are invalidated.
- Neighbor rebuild: compiled neighbor arrays are refreshed.
- Output/checkpoint: may force interaction evaluation.
