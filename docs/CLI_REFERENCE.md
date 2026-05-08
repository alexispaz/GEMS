# GEMS CLI Reference

Related documents: [input language](INPUT_LANGUAGE.md), [architecture](ARCHITECTURE.md), and [developer guide](DEVELOPER_GUIDE.md).

Notation: `g` is a saved group index, `vec3` is three numbers, and `:label` is optional unless stated.

## Core Commands

| Command | Syntax | Effect |
| --- | --- | --- |
| `license` | `license` | Print license. |
| `warranty` | `warranty` | Print warranty. |
| `if` | `if int command` | Execute one command token if `int == 1`. |
| `exit`, `stop` | `exit` | End current command-processing loop. |
| `cycle` | `cycle` | Cycle signal, also recognized in blocks. |
| `dimension` | `dimension int` | Legacy/unclear dimension metadata. |
| `temp` | `temp T` | Set global temperature for `kt` conversion. |
| `walltime` | `walltime` | Print elapsed CPU time. |
| `help` | `help topic` | Stubbed. |

## Create / Select / Destroy

Creation subcommands:

| Syntax | Effect |
| --- | --- |
| `+ atom vec3` | Create atom in `sys`. |
| `>+ atom vec3` | Create atom and replace selection. |
| `^+ atom vec3` | Create atom and add to selection. |
| `... read file [frame]` | Read atoms from file. |
| `... reply vec3 n` | Replicate current selection. |
| `... fill n radius [min vec3 max vec3]` | Fill atoms with spacing. |
| `... fillpbc n radius [min vec3 max vec3]` | Fill with PBC option. |
| `... fillh` | Add hydrogens to selected carbons. |

Selection subcommands:

| Syntax | Effect |
| --- | --- |
| `> all` | Select all system atoms. |
| `> group g` | Select saved group `g`. |
| `> element E` | Select element. |
| `> atom i [j]` | Select atom/range. |
| `^> sphere r [center]` | Intersect current selection with sphere. |
| `^> cylinder r axis [center]` | Intersect with cylinder. |
| `^> xrange/yrange/zrange lo hi` | Coordinate band. |
| `^> rectangle min max` | Rectangular region. |
| `^> above p1 p2 p3`, `^> below p1 p2 p3` | Plane filters. |
| `^> near vec3`, `^> far vec3` | Point-distance selectors. |
| `^> bo_range group lo hi [r1 r2]` | Bond-order range. |
| `^> bo_min group [r1 r2]`, `^> bo_max group [r1 r2]` | Bond-order extrema. |
| `^> random n` | Random subset. |

Destroy:

| Syntax | Effect |
| --- | --- |
| `- selector` | Destroy atoms selected from `sys`. |
| `^- selector` | Destroy atoms selected from current selection. |
| `^~ selector` | Remove selected atoms from current selection without destroying them. |

## Groups

| Syntax | Effect |
| --- | --- |
| `group g add` | Attach current selection to saved group `g`. |
| `group g clean` | Clear saved group `g`. |
| `clean creation` | Select `sys` and try to destroy atoms no longer referenced. |

## Set

All `set` subcommands operate on current selection.

| Syntax | Effect |
| --- | --- |
| `set mass m` | Set mass. |
| `set cm_vel vec3` | Set center-of-mass velocity. |
| `set pbc bool bool bool` | Set atom PBC and allocate/deallocate ghosts. |
| `set tempgdist T` | Assign Gaussian velocity distribution. |
| `set cm_pos vec3`, `set cg_pos vec3` | Move group center. |
| `set sigma s` | Set sigma. |
| `set minpos axis value`, `set maxpos axis value` | Shift bounding coordinate. |
| `set posdist file [frame]`, `set veldist file [frame]` | Read positions/velocities. |
| `set element symbol` | Set element. |
| `set min_interdist r` | Enforce minimum distance. |
| `set align`, `set align_mass`, `set alignxy_mass` | Alignment helpers. |
| `set align_by g`, `set align_massxy_by g` | Align relative to group. |
| `set align_tutor vec3 vec3` | Tutor alignment. |
| `set rotateaxis angle p v` | Rodrigues rotation, angle in degrees. |
| `set rotate x/y/z angle [origin]` | Axis rotation. |
| `set expand vec3 [origin]` | Scale positions. |
| `set move vec3` | Translate positions. |
| `set do_pbc` | Fold selected atoms through PBC. |

## Box

| Syntax | Effect |
| --- | --- |
| `box size x y z` | Set box. |
| `box expand ax ay az` | Scale box. |
| `box move` | Remove total CM velocity. |
| `box mic bool` | `T`: minimum image; `F`: ghost mode. |
| `box tsize` plus matrix lines | Non-cubic transform path; limited by source comments. |

## Interactions

```gms
interact [:label] g kind subkind params...
interact [:label] g < h kind subkind params...
```

Kinds:

| Kind | Subkinds / Parameters |
| --- | --- |
| `pair` | `lj eps sigma rcut`, `plj eps sigma`, `slj eps sigma rcut`, `wca eps sigma`, `rwca eps sigma start`, `cuw eps sigma start radius coordinate rcut`, `sho k r0`, `shocm k r0`, `sm1 p1 p2 energy`. |
| `field` | `halfsho_plane`, `sho2d`, `voter2d`, `oscar2d`, `leiva1d`, `pozoa1d`, `sho_plane`, `sho_line`, `lucas1d`; parameters are delegated to `field_cli`. |
| `bias` | `mcammon`, `compress_below`, `compress`. |
| `graph` | `subgraphs`, `blocks`. |
| `tb` | `read file` or `a eps p q r0 rci rce`. |

Examples:

```gms
interact 1 pair lj 1.0 0.89 2.0
interact 2 < 3 pair lj 0.1 3.4 200
interact 1 tb read test.prm
```

## Evolve

```gms
evolve [:label] algorithm params...
```

| Algorithm | Parameters |
| --- | --- |
| `v_verlet` | none |
| `pred_corr_5` | none |
| `brownian` | two real parameters |
| `ermak` | two real parameters; source reads second then first into setup call |
| `ermak_x`, `ermak_y`, `ermak_z` | directional Ermak |
| `gcmc` | five real/integer parameters as read by `set_gcmc` |
| `piston_lkd` | five parameters |
| `piston_lgf` | axis plus five parameters |
| `secfile file` | Read positions/velocities from file |
| `scalvel value each` | Velocity scaling |
| `scalvel_after value each after` | Delayed scaling |
| `andersen freq [window]` | Andersen thermostat |
| `voter bool` | Hyperdynamics time scaling; requires Ermak |
| `clean` | Destroy integration objects |

## Output

Immediate:

| Syntax | Effect |
| --- | --- |
| `out posxyz [file]` | Write current selection. |
| `out state` | Run `interact(.true.)` and write force output. |
| `out float fmt`, `out integer fmt` | Set output formats. |

Persistent:

| Syntax | Effect |
| --- | --- |
| `outfile [:label] name file` | Set output filename. |
| `outfile [:label] each n` | Sample interval. |
| `outfile [:label] flush on/off` | Flush control. |
| `outfile [:label] prom n` | Average over `n`. |
| `outfile [:label] ddda n` | DDDA error analysis. |
| `outfile [:label] at common/eparts/forcefield/hd` | Write phase. |
| `outfile [:label] cols col [group] ...` | Numeric columns. |
| `outfile [:label] pos|poscr|vel|fce|pes|pose|dist|charge group` | Structured output. |
| `outfile [:label] calc` | Calculation output. |

Common columns include `time`, `step`, `box`, `eparts`, `bias`, `energy`, `epot`, `ecin`, `temp`, `pressure`, `virial`, `girrad`, `cvs`, `inercia`, `covariance`, `mainaxis`, and `ptriaxial`.

## Print / Get

`print query` writes a value. `getin var query` stores the same result.

Queries include `dm_steps`, `selnat`, `temp`, `cm_vel`, `cm_pos`, `cm_diff`, `cm_diff2`, `rg_pos`, `minpos`, `maxpos`, `ptriaxial`, `mayordist`, `border`, `index`, `rang`, `ranu`, `max`, `min`, `int`, `floor`, `ceil`, and `box`.

## Checkpoint, PRNG, MPI, Time

| Command | Syntax |
| --- | --- |
| `checkpoint` | `write file`, `read file`, `on`, `off`, `each n`. |
| `prng` | `lcg`, `std`, `seed n`, SPRNG-only `seed_first`, `stream`. |
| `mpi` | `like pc tpc` in serial builds; `only pc command`. |
| `time` | `cero`, `step dt`. |
| `element` | `add symbol mass`. |

## Specialized Workflows

NEB/minimization:

- `neb_clean`, `neb_nimg n`, `neb_img i`, `neb_run iterations k`
- `minvol f1 f2 [bool]`
- `lbfgs [bool]`

Metadynamics/hyperdynamics/parallel tempering:

- `setwtmd1d`, `setpos1d`, `setwall1d`, `wtdcm`, `wtx`, `wtpos1d`
- `setwtmd2d`, `setwall2d`, `setwallauxcore`, `setwallauxshell`
- `wtdcmrgco`, `wtxy`, `wtdcmrgau`, `wtdcmrgtotal`, `dmcv`
- `hybrid_wb`, `hybrid_ea`
- MPI-only `partemp`, `partemp2`

Use `examples/WTMD`, `examples/hyperdynamics`, `examples/RepExchange`, `examples/NEB`, and `examples/optimize` for working argument examples.
