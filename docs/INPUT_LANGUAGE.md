# GEMS Input Language

Related documents: [CLI reference](CLI_REFERENCE.md), [architecture](ARCHITECTURE.md), and [dataflow](DATAFLOW.md).

## Overview

GEMS input is a line-oriented scripting language. Each non-comment line is parsed into items, the first item is interpreted as a command, and the remaining items are consumed by that command or its subcommand.

Typical file structure:

```gms
box size 5.78 5.78 5.78
>+ atom 1 1 1
set element H
> all
group 1 add
interact 1 pair lj 1.0 0.89 2.0
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
evolve v_verlet
dinamica 100
```

## Parser Flow

`Input_Parsing.F90::read_line` performs:

1. physical line read, with include/revert support;
2. line continuation handling;
3. tab replacement and leading-space trim;
4. directive handling for `#include`, `#concat`, `#revert`, `#width`, `#echo`;
5. comment skipping;
6. variable assignment handling;
7. variable and math expansion;
8. item tokenization;
9. block handling for `bloque ...`.

`Main.F90` then reads the first item with `readl` and calls `execute_command`.

## Items And Comments

- Items are separated by spaces, tabs, or commas.
- Single quotes keep text as one item.
- Parentheses are stripped as comments, including nested parenthetical comments.
- Lines beginning with `# ` are comments.
- Lines beginning with `#include`, `#revert`, etc. are directives.

Example:

```gms
set move 1 0 0  (translate current selection)
outfile :f1 name 'Energy with spaces.dat'
```

## Variables And Math

Assignments:

```gms
name=value
name:={sqrt(0.5)}
```

- `=` stores text.
- `:=` expands variables and evaluates math.
- Variables are referenced as `$name` or `$name$`.
- `$rnd` expands to a random number.
- User variable names cannot start with `_`.

Math expressions are evaluated by the bundled `fparser` library.

Examples:

```gms
a0:=1
d0:={sqrt(0.5)}
box size $a0$ $a0$ $a0$
>+ atom {0.5*$a0$} 0 0
```

Startup variables include `jobname`, `i`, `ci`, `boxx`, `boxy`, `boxz`, `time`, `step`, `dt`, and read-only `pi`. MPI builds also define MPI rank variables.

## Blocks

Blocks begin with `bloque` and end with `fin`.

| Syntax | Effect |
| --- | --- |
| `bloque repeat n` | Execute following block for `1..n`. |
| `bloque repeat a b [step]` | Execute range; reverse if `a>b`. |
| `bloque save k` | Store a block in memory slot `k`. |
| `bloque load k` | Execute saved block once. |
| `bloque load_each each k` | During dynamics, execute saved block `k` every `each` steps. |
| `bloque load_off` | Disable scheduled block execution. |
| `bloque load_silent bool` | Set load-silence flag. |

During repeat, `$i` and `$ci` are updated. `cycle` and `exit` are handled specially inside block execution.

## Selection And Creation Model

The current CLI creates atoms directly in `sys` and manages a current selection `gsel`.

| Command | Meaning |
| --- | --- |
| `+ ...` | Create atoms into `sys`. |
| `>+ ...` | Create atoms and replace current selection with them. |
| `^+ ...` | Create atoms and union them into current selection. |
| `> ...` | Select from `sys`, replacing current selection. |
| `^ ...` | Select from `sys`, union into current selection. |
| `^> ...` | Select from current selection, replacing it. |
| `^~ ...` | Remove selected atoms from current selection. |
| `- ...` | Destroy atoms selected from `sys`. |
| `^- ...` | Destroy atoms selected from current selection. |

Examples:

```gms
>+ read coords.xyz
> all
> element Au
^> sphere 6.0 0 0 0
^~ bo_range -1 0.9 1.0 1.7 1.9
```

Saved memory groups use:

```gms
group 1 add
group 1 clean
```

## Interactions

Interactions are `ngroup` objects. Syntax:

```gms
interact [:label] group kind subkind params...
interact [:label] group < other_group kind subkind params...
```

Examples:

```gms
interact 1 pair lj 1.0 0.89 2.0
interact 2 < 3 pair lj 0.1 3.4 200
interact :gr 1 graph subgraphs 3.0
```

The first group is `ref`; the optional group after `<` is the neighbor candidate group `b`. If omitted, `b=ref`.

Interaction setup immediately calls `test_update`, so declaring interactions can build/update ghost, cell, and neighbor structures.

## Side Effects To Know

- `out state` and checkpoint read/write call `interact(.true.)`.
- `interact` calls `test_update`, which can update PBC, ghosts, cells, and neighbor lists.
- `set pbc` allocates or deallocates atom ghost arrays.
- `box mic F` sets `useghost=.true.`; `box mic T` uses minimum-image mode.
- Integration steps invalidate all group caches.
- Group attach/detach invalidates group properties and can invalidate indexes in `igroup` descendants.

## Inconsistencies / Ambiguities

- `help` is currently stubbed.
- `dimension` still appears to write dimension metadata but the real dimension is compile-time controlled by `dm`.
- Command names mix symbols, English, and Spanish (`dinamica`, `cero`).
- `table` and `mtable` have mostly commented-out subcommands.
- Math capabilities depend on `fparser`; this document does not enumerate every function.
