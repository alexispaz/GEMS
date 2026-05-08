# GEMS Architecture Notes

Status: fresh analysis of `main` at `af392e6` on branch `codex-docs`.

Related documents: [memory model](MEMORY_MODEL.md), [dataflow](DATAFLOW.md), [hot paths](HOTPATHS.md), [input language](INPUT_LANGUAGE.md), [CLI reference](CLI_REFERENCE.md), and [developer guide](DEVELOPER_GUIDE.md).

## Overview

GEMS is a command-driven molecular simulator written in Fortran. The current code is organized around dynamic groups of atom objects, indexed group variants for faster lookup, neighbor groups for force kernels, and a line-oriented input interpreter.

The current architecture more explicitly separates these layers than older versions:

```text
input commands
    -> editable groups and atom objects
    -> indexed groups / cell groups / neighbor groups
    -> compute-heavy interaction, integration, output, and analysis routines
```

The conceptual model is partly implemented:

```text
editable dynamic representation
    -> compiled dense representation
    -> compute-heavy kernels
```

The dense representation is not a single global compiled state. It is distributed across `igroup`, `cgroup`, and `ngroup`.

## Main Modules

| Module | Role |
| --- | --- |
| `Main.F90` | Program entry point, startup, global registries, input loop. |
| `Program_Types.F90` | Box, time, PBC-distance helpers, global run variables. |
| `Groups.F90` | Core `atom`, `group`, `igroup`, system group, ghost management, group property invalidation. |
| `Cells.F90` | `cgroup`, linked-cell tessellation and sorting for indexed groups. |
| `Neighbor.F90` | `ngroup`, neighbor-list storage, neighbor update policy, neighbor search dispatch. |
| `Interaction.f90` | Interaction creation and per-step force/energy dispatch. |
| `Integration.F90` | Integrator objects and integration step dispatch. |
| `Programs.f90` | High-level run loops such as `dinamic`. |
| `CLInterpreter.F90` | Top-level input commands and command-specific subinterpreters. |
| `Input_Parsing.F90` | Line reading, directives, variables, math expansion, tokenization, blocks. |
| `Pairs.F90`, `Fields.f90`, `TB.F90`, `Bias.f90`, `Graphs.F90`, `ForceField.F90` | Concrete interaction and analysis implementations. |
| `Calc.F90` | Calculation objects such as Widom-style calculations and output support. |
| `Output.F90`, `Checkpoint.F90`, `CVS.F90`, `NEB.f90`, `Metadynamics.f90`, `Hyperdinamics.f90` | Output, restart, CVs, and specialized workflows. |

## Core Structure Types

### `atom`

Defined in `Groups.F90`. It owns mechanical vector pointers (`pos`, `vel`, `force`, `acel`), element/scalar fields, PBC flags, ghost links, and group membership metadata.

An atom records group membership through sorted `gr(:)` group ids and matching `id(:)` per-indexed-group ids.

### `group`

Editable atom membership represented by a linked list `alist`. It is used for selections, saved groups, and general dynamic membership.

Groups are registered in `gindex`, have stable `id` values, and expose invalidation methods for cached group properties.

### `igroup`

Extends `group` with a compact-ish array of atom pointers:

```fortran
type(atom_ap), allocatable :: a(:)
integer :: amax
```

`igroup` gives atoms a group-local id and supports cleanup/compaction after detaches. This is the main bridge from linked membership to indexed access.

### `cgroup`

Extends `igroup` with linked-cell storage:

```fortran
integer, allocatable :: head(:,:,:)
integer, allocatable :: next(:)
```

It tessellates the simulation box and sorts indexed atoms into cells for neighbor search.

### `ngroup`

Extends `igroup` and represents an interaction group with neighbors. It owns:

- `ref`: atoms for which forces are computed.
- `b`: candidate neighbor group, as a `cgroup`.
- `nn(:)` and `list(:,:)`: dense neighbor-list arrays.
- procedure pointers for neighbor search strategy.
- deferred `interact` and `cli` procedures implemented by concrete interaction modules.

`ngroup` instances are registered in `ngindex` and executed in order.

## Execution Flow

1. `Main.F90` initializes RNG, elements, variables, groups, registries, and parser options.
2. The input loop calls `read_line`, reads the first token, then calls `execute_command`.
3. Commands create atoms directly into `sys`, modify `gsel`, save `gr(i)`, configure interactions/output/integrators, and run workflows.
4. `dinamic` performs an initial `interact(.false.)`, then repeats integration A, interaction, integration B, output/checkpoint/block hooks.
5. `interact` calls `test_update`, zeroes forces/energies for system and ghosts, dispatches every `ngroup`, dispatches bonded force-field groups, writes optional output, and evaluates CV force aggregation.

## Current Dynamic-To-Dense Pipeline

Implemented pieces:

- `group%alist`: editable arbitrary selections.
- `igroup%a(:)`: indexed atom pointer array compiled from group membership.
- `cgroup%head/next`: linked-cell compact search structure.
- `ngroup%nn/list`: dense neighbor-list storage consumed by interaction kernels.
- `group_switch_vectorial`: optional contiguous group vectors for selected vector algorithms.

Missing or incomplete pieces:

- There is no single immutable compiled simulation state.
- Atom payloads remain pointer-rich even when neighbor lists are dense.
- Ghost activation/deactivation mutates group membership at runtime.
- Invalidation is distributed through group flags, `listed/tessellated` flags, and displacement checks.

## Architectural Invariants

- All groups must be indexed in `gindex` and have an id.
- Atom group membership is mirrored: group linked-list membership and atom `gr/id` membership must stay consistent.
- `igroup` indices are group-local; `atom%gid(g)` maps an atom to its id inside a specific indexed group.
- `ngroup%list(i,m)` stores neighbor ids in the `ngroup` index space, not raw system atom ids.
- `ngroup%ref` and `ngroup%b` must be attached before `ngroup` itself is used for neighbor sorting.

## Open Questions

- The intended long-term split between `group`, `igroup`, `cgroup`, and `ngroup` is clearer than before, but still not documented in source-level invariants.
- `group_switch_vectorial` is marked TODO/delete in comments but still present and used by algorithms; its future role is uncertain.
- Interaction creation no longer has the older explicit `nomoreigs` freeze. Runtime behavior now allows `interact` declarations to call `test_update` immediately; whether late interaction creation is intended needs project-owner confirmation.
