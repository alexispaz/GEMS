# Hot Paths And HPC Readiness

Related documents: [architecture](ARCHITECTURE.md), [memory model](MEMORY_MODEL.md), [dataflow](DATAFLOW.md), and [developer guide](DEVELOPER_GUIDE.md).

No profiling was run for this pass. This is a static hot-path inventory.

## Likely Hot Paths

### `interact`

File: `src/Interaction.f90`

- Called before and during dynamics.
- Always calls `test_update`.
- Zeros forces/energies over `sys` and `ghost`.
- Dispatches every `ngroup` polymorphically.
- Accumulates energies and CV hooks.

### Neighbor Update

Files: `src/Neighbor.F90`, `src/Cells.F90`

- `test_update` runs each interaction call.
- `cgroup%tessellate` and `cgroup%sort` prepare linked cells.
- `ngroup_verlet` and `ngroup_cells` build dense neighbor arrays.
- OpenMP tasks are used inside neighbor builders.

### Pair/TB/Field Kernels

Files: `src/Pairs.F90`, `src/TB.F90`, `src/Fields.f90`, `src/Bias.f90`, `src/Graphs.F90`

- Pair/TB kernels consume `ngroup%nn/list`.
- Field-like kernels may iterate directly over linked group membership.
- TB owns additional dense parameter/band arrays and can be heavier per neighbor.

### Group Property Queries

Files: `src/Inq_Properties.F90`, `src/Set_Properties.f90`, `src/Output.F90`

- Many property queries traverse group linked lists.
- Output frequency can make these paths important.

## Allocation Hot Spots

Potentially frequent:

- Atom allocation in create/read/fill workflows.
- `set_pbc` ghost array allocation/deallocation.
- ghost attach/detach in `pbchalfghost` / `ghost_switch`.
- `igroup` growth via `move_alloc`.
- `cgroup` head/next allocation when tessellation changes.
- `ngroup` neighbor-list growth when attaching atoms.

Mostly setup/workflow-specific:

- interaction object allocation.
- TB parameter matrices.
- CV, NEB, output, DDDA arrays.

## Locality Risks

Current force-kernel access generally follows:

```text
ngroup%list(i,m)
    -> ngroup%a(j)%o
    -> atom%pos / atom%force pointer targets
```

Dense neighbor arrays help, but atom payload remains pointer-rich. This limits cache locality, compiler vectorization, and future GPU offload.

## GPU/OpenMP Preparation

Already useful:

- `igroup` and `ngroup` define stable local ids.
- neighbor lists are dense arrays.
- cell structures are separate from interaction formulas.
- displacement-based rebuild policy exists.

Needs systematization:

- explicit compiled atom storage policy;
- explicit ghost/halo arrays;
- kernel contracts that separate read-only positions from force reductions;
- removal of linked-list traversal from compute kernels;
- controlled lifecycle/dirty flags for compiled state.

## Documentation-Only Opportunities

Before optimizing:

- inventory each concrete `ngroup` by required atom fields;
- document `ngroup%list` invariants and index ranges;
- document ghost ownership and attachment rules;
- mark commands that force neighbor rebuilds or interaction evaluation;
- add developer notes around `group_switch_vectorial`.
