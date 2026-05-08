# Developer Guide

Related documents: [architecture](ARCHITECTURE.md), [memory model](MEMORY_MODEL.md), [dataflow](DATAFLOW.md), [hot paths](HOTPATHS.md), [input language](INPUT_LANGUAGE.md), and [CLI reference](CLI_REFERENCE.md).

## Build Systems

The current README recommends Meson, with Autotools still available.

Meson:

```sh
meson setup build
meson install -C build
```

Autotools:

```sh
autoreconf -fi
./configure
make
make install
```

The example test harness expects an installed executable at `usr/bin/gems`:

```sh
cd examples
./test.sh
```

## Documentation Preview

The documentation site is configured with MkDocs and Material for MkDocs:

```sh
python -m venv .venv-docs
. .venv-docs/bin/activate
python -m pip install -r requirements-docs.txt
mkdocs serve
```

The MkDocs configuration is intentionally local-preview only. It does not configure deployment or GitHub Pages.

## Source Orientation

- Parser and CLI: `Input_Parsing.F90`, `CLInterpreter.F90`.
- Mutable state: `Groups.F90`, `Program_Types.F90`.
- Indexed/cell/neighbor infrastructure: `Groups.F90`, `Cells.F90`, `Neighbor.F90`.
- Runtime loops: `Main.F90`, `Programs.f90`, `Interaction.f90`, `Integration.F90`.
- Concrete interactions: `Pairs.F90`, `Fields.f90`, `TB.F90`, `Bias.f90`, `Graphs.F90`, `ForceField.F90`.
- Analysis/output/workflows: `Output.F90`, `Calc.F90`, `CVS.F90`, `NEB.f90`, `Metadynamics.f90`, `Hyperdinamics.f90`, `DDDA.F90`.

## Safe Change Rules

- Keep physical behavior unchanged unless explicitly validating physics.
- Preserve CLI syntax or add compatibility aliases.
- Update `docs/CLI_REFERENCE.md` when command behavior changes.
- Update examples/tests when adding user-visible features.
- Avoid changing ownership of `atom%pos/vel/force/acel` without auditing `group_switch_vectorial`.

## Adding A CLI Command

1. Add a `case` in `execute_command` or an existing `*_commands` routine.
2. Use existing readers: `readl`, `reada`, `readi`, `readf`, `readb`.
3. If it mutates atoms/groups, trigger appropriate invalidation.
4. Add a working example when possible.
5. Document syntax and side effects.

## Adding An Interaction

1. Define a concrete type extending `ngroup`.
2. Implement `interact` and `cli`.
3. Add allocation in the relevant factory (`pair_new`, `field_new`, `bias_new`, `graph_new`, or `interact_new`).
4. Attach `ref`, `b`, and the `ngroup` consistently.
5. Set `rcut` via `setrc` if neighbor lists are required.
6. Document parameters in `CLI_REFERENCE.md`.

## Memory Work To Clarify Before Refactors

- Formalize whether `group_switch_vectorial` remains supported.
- Document `ngroup%list` index invariants in source comments.
- Decide whether late interaction creation after dynamics has begun is intended.
- Define future dense atom storage as a separate compiled view before GPU work.

## Logging Documentation Work

Record substantial documentation/analysis work in `docs/dev/AGENT_LOG.md` with inspected files, decisions, uncertainty, and tests.
