# GEMS Memory Model

Status: fresh analysis of current `main` at `af392e6`.

Related documents: [architecture](ARCHITECTURE.md), [dataflow](DATAFLOW.md), [hot paths](HOTPATHS.md), and [developer guide](DEVELOPER_GUIDE.md).

## Conceptual Model

The current code supports the intended model:

```text
editable dynamic representation
    -> compiled dense representation
    -> compute-heavy kernels
```

but implements it as a set of layered structures rather than a monolithic compiled-state object.

Editable state:

- `group%alist`: linked-list membership.
- `gsel`, `gr(:)`, `sys`, `ghost`, and workflow-specific groups.
- Per-atom fields and PBC/ghost membership.

Compiled/indexed state:

- `igroup%a(:)`: dense-ish atom pointer index.
- atom `gr(:)`/`id(:)` membership arrays.
- `cgroup%head(:,:,:)` and `next(:)`.
- `ngroup%nn(:)` and `list(:,:)`.

Compute-heavy state:

- concrete `ngroup` interaction objects (`lj`, `smatb`, field/bias/graph types).
- integration groups.
- output accumulators and analysis/CV arrays.

## Ownership

### Atom Objects

Atoms are heap-allocated and attached to groups. The primary owner is not an isolated atom container; lifetime is controlled by group membership.

- `atom%dest` detaches the atom from every group before deallocating internal arrays.
- `group%destroy_all` calls `a%dest()` and deallocates atom objects.
- `group%try_destroy_all` detaches atoms and deallocates only atoms with no remaining group references.
- Ghost atoms are either preallocated in `atom%ghost(1:7)` when PBC is enabled or allocated by `new_ghost` in the older full ghost path.

### Groups

`group` owns its membership list nodes, not necessarily the atom objects. Attaching an atom also registers the group id in the atom. Detaching removes both sides.

`gindex` stores pointers to groups. `group_destroy` removes the group from `gindex` but does not renumber group ids.

### Indexed Groups

`igroup` owns `a(:)` and `limbo`.

- `a(:)` stores atom pointers, not atom payload.
- `amax` is the active scanned extent.
- Detach can null slots or put atoms in `limbo`.
- `igroup_update_index` compacts the pointer index and rewrites atom-local ids.

### Cell Groups

`cgroup` owns cell-list arrays:

- `head(:,:,:)`: first atom index per cell.
- `next(:)`: next atom index in cell chain.

It does not own atoms; it sorts existing indexed atoms.

### Neighbor Groups

`ngroup` owns:

- internal `ref` group.
- internal candidate `b` `cgroup`.
- inherited `igroup` index.
- neighbor arrays `nn(:)` and `list(:,:)`.

Concrete interaction objects extend `ngroup` and may own additional arrays, for example TB parameter matrices and band arrays.

## Atom Field Storage

`atom%pos`, `vel`, `force`, and `acel` are Fortran pointers. By default, `atom_allocate` allocates each as a separate `dm`-length heap target.

`group_switch_vectorial` can replace per-atom vector allocations with slices into group-owned contiguous arrays:

```text
atom vector allocations
    -> group%pp/pv/pf/pa contiguous arrays
    -> atom vector pointers redirected to slices
```

This improves contiguity for a chosen group but changes ownership of the atom vector targets. `group_switch_objeto` reverses it. The source comments mark this area as a workaround/TODO, so future refactors should treat it as fragile.

## Invalidation And Recompile

Property invalidation:

- `gindex_all_changed`, `gindex_pos_changed`, `gindex_vel_changed`, and `gindex_epot_changed` broadcast cache invalidation to all registered groups.
- `group_attach_atom` and `group_detach_link` call `all_changed`.
- integration steps call `gindex_all_changed`.

Indexed group invalidation:

- `igroup_detach_atom` may null slots and trigger `igroup_update_index` when holes exceed `aupd`.
- `igroup%update` signals extended types that internal arrays may need rebuild.

Cell-list invalidation:

- `cgroup%tessellated` tracks whether cell storage is valid.
- `cgroup_attach_atom` tries to sort new atoms incrementally when tessellated; failures mark it stale.
- `cgroup_detach_atom` can unsort or mark stale.
- `test_update` tessellates and sorts `ngroup%b` before checking neighbor rebuilds.

Neighbor-list invalidation:

- `ngroup%listed` indicates whether `nn/list` are valid.
- `test_update` rebuilds if a neighbor group is not listed.
- It also rebuilds when the sum of the two largest displacements exceeds `nb_dcut`.
- `update` records `pos_old` and calls each ngroup's selected neighbor builder.

Ghost invalidation:

- `set_pbc` allocates/deallocates ghost arrays.
- `box mic F` enables `useghost`; `test_update` calls `pbchalfghost(maxrcut)`.
- `ghost_switch` attaches/detaches active ghost atoms from every ghost-enabled group that contains the prime atom.

## Dynamic Vs Compute-Oriented Structures

Dynamic/flexible:

- `group`, `gsel`, `gr(:)`.
- atom group membership arrays.
- ghost activation/deactivation.
- parser variables and labels.

Compute-oriented:

- `igroup%a(:)` for indexed atom access.
- `cgroup` cell arrays.
- `ngroup%nn/list`.
- concrete interaction arrays.
- output/CV/NEB/TB work arrays.

Partially compute-oriented but still pointer-rich:

- force kernels use neighbor ids but dereference atom pointers and atom vector pointers.
- OpenMP tasks exist in neighbor search, but atom payload layout is not structure-of-arrays.

## Locality / Cache Implications

Good locality:

- neighbor counts and lists are dense integer arrays.
- cell heads and next arrays are dense.
- group vectorization can make selected vector fields contiguous.

Poor or uncertain locality:

- atom objects are individually allocated.
- atom mechanical fields are separate pointer targets unless vectorized.
- linked-list traversal remains common in setup, properties, output, and some kernels.
- `ngroup%list -> ngroup%a -> atom -> pos/force` still has multiple pointer hops.

## OpenMP/GPU Readiness

Promising pieces:

- explicit indexing layers (`igroup`, `ngroup`).
- dense neighbor-list arrays.
- clear separation between neighbor update and interaction dispatch.
- linked-cell data structure.

Friction points:

- object/pointer-heavy atom payload.
- polymorphic interaction dispatch.
- ghost membership mutates groups at runtime.
- global module state and broad invalidation.
- no explicit device-ready compiled state or storage policy.

Before GPU work, document and stabilize a kernel contract: dense input arrays, force/energy outputs, ghost/halo ownership, and rebuild policy.
