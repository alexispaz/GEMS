# Agent Log

## 2026-05-08 09:39:42 -0300

### Scope

Fresh documentation pass after synchronizing local `main` to `origin/main` at `af392e6` and creating branch `codex-docs`.

### Files Inspected

- `README.md`
- `meson.build`
- `src/meson.build`
- `src/Makefile.am`
- `src/Main.F90`
- `src/Program_Types.F90`
- `src/Groups.F90`
- `src/Cells.F90`
- `src/Neighbor.F90`
- `src/Interaction.f90`
- `src/Programs.f90`
- `src/Integration.F90`
- `src/CLInterpreter.F90`
- `src/Input_Parsing.F90`
- `src/Pairs.F90`
- `src/Fields.f90`
- `src/TB.F90`
- `src/Bias.f90`
- `src/Graphs.F90`
- `src/Calc.F90`
- `examples/test.sh`
- representative `.gms` files under `examples/`

### Documentation Created

- `docs/ARCHITECTURE.md`
- `docs/MEMORY_MODEL.md`
- `docs/DATAFLOW.md`
- `docs/HOTPATHS.md`
- `docs/INPUT_LANGUAGE.md`
- `docs/CLI_REFERENCE.md`
- `docs/DEVELOPER_GUIDE.md`

### Important Findings

- Current `main` has a different architecture from the obsolete previous pass.
- `Program_Types.F90` now mostly owns box/time/PBC helpers, while `Groups.F90` owns `atom`, `group`, `igroup`, `sys`, and ghost management.
- Dense/indexed state is layered: `igroup%a`, `cgroup%head/next`, and `ngroup%nn/list`.
- `Neighbor.F90` now defines `ngroup` rather than the previous `intergroup` design.
- The CLI now uses `+`, `>+`, `^+`, `>`, `^`, `^>`, `^~`, `-`, and `^-` for create/select/destroy flows.
- `interact` syntax is now `interact [:label] group [< group] kind subkind ...`.
- Meson is now the recommended build system; Autotools remains available.

### Decisions Taken

- Did not reuse previous conclusions automatically.
- Kept existing `docs/ddda.md` and `docs/types.md` intact.
- Documented ambiguous behavior explicitly.
- Did not modify Fortran source.

### Uncertainties

- Whether late interaction creation after dynamics has begun is intended.
- Long-term support status of `group_switch_vectorial`.
- Full semantics of specialized metadynamics/hyperdynamics commands beyond source comments and examples.
- Complete function set exposed by bundled `fparser`.

### Tests Executed

- Not run yet; documentation-only changes so far.

## 2026-05-08 10:01:34 -0300

### Scope

MkDocs local preview configuration.

### Files Inspected

- `docs/ARCHITECTURE.md`
- `docs/MEMORY_MODEL.md`
- `docs/DATAFLOW.md`
- `docs/HOTPATHS.md`
- `docs/INPUT_LANGUAGE.md`
- `docs/CLI_REFERENCE.md`
- `docs/DEVELOPER_GUIDE.md`
- `docs/dev/AGENT_LOG.md`

### Configuration Added

- Added root `mkdocs.yml`.
- Set `docs_dir: docs`.
- Set Material for MkDocs with:
  - `theme.name: material`

### Navigation Decision

The MkDocs navigation lists only the main user/developer documentation pages:

1. Architecture
2. Memory model
3. Data flow
4. Hot paths
5. Input language
6. CLI reference
7. Developer guide

Existing auxiliary documents such as `docs/ddda.md`, `docs/types.md`, and `docs/dev/AGENT_LOG.md` remain in `docs/` but are intentionally omitted from the primary navigation for now.

### Deployment

- No deployment configured.
- No GitHub Pages configuration added.

### Tests Executed

- `mkdocs --version` succeeded with MkDocs 1.6.1.
- `mkdocs build --strict --site-dir /tmp/gems-mkdocs-site` failed because the local environment has MkDocs installed but not the `material` theme.
- Added `requirements-docs.txt` with `mkdocs-material` so local preview can be installed in a project-local virtual environment instead of globally.
- Created a temporary virtual environment at `/tmp/gems-docs-venv` for validation only.
- Installed `mkdocs-material` into the temporary virtual environment after sandboxed network access failed.
- `/tmp/gems-docs-venv/bin/mkdocs build --strict --site-dir /tmp/gems-mkdocs-site` succeeded.
- MkDocs reported that `ddda.md`, `types.md`, and `dev/AGENT_LOG.md` are not in `nav`; this matches the documented navigation decision.
