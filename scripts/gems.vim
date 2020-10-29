" Vim syntax file
" Language:	GMD
" Maintainer:	alexis paz <apaz@fcq.unc.edu.ar>
" Extensions:   *.gmd
" Comment:      GMD is a input script language for General Molecular Dynamic program.

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

" Ignore case
syn case ignore

" A bunch of useful gmd keywords

"syn match   gmdIdentifier

syn keyword gmdBlock       bloque 
syn keyword gmdBlock       fin
syn keyword gmdBlock       repeat

syn keyword gmdStop        stop

syn keyword gmdIdentifier  out outfile box neigh time
syn keyword gmdSelected    set interact under evolve
syn match   gmdAdd         "\p* add"
syn match   gmdSelected    "\s<\s"
syn match   gmdBCrea       "{\p\{-}}"
syn match   gmdSelected    "\$\p\{-}\$"

"syn keyword gmdGsel        sys group 
syn match   gmdGsel        "group\b*\d"
syn match   gmdGsel        ">"
syn match   gmdGsel        "+"
syn match   gmdGsel        ">\s\p*"
syn match   gmdGsel        "+\s\p*"
syn match   gmdGsel        ">>\s\p*"
syn match   gmdGsel        "-\s\p*"
syn match   gmdBCrea       "<\s\p*"
syn match   gmdBCrea       "^\."

" \zs and \ze set the start and the final of the submatch to higlight
syn match   gmdSelected    "\d\s\s*\zs<>\ze\s"
syn match   gmdSelected    "\d\s\s*\zs>\ze\s"

syn keyword gmdStatement   partemp dinamica hd lbfgs cg lp

"syn match gmdBoolean	"\ \s*\(T\|F\)\s* "
"
"
"syn match   gmdStatement   "flag-\(non45\|acuteangle\|offgrid\)"
"syn match   gmdStatement   "text-pri-only"
"syn match   gmdStatement   "[=&]"
"syn match   gmdStatement   "\[[^,]*\]"
"syn match   gmdstatement   "^ *\(sel\|width\|ext\|enc\|area\|shrink\|grow\|length\)"
"syn match   gmdstatement   "^ *\(or\|not\|and\|select\|size\|connect\|sconnect\|int\)"
"syn match   gmdstatement   "^ *\(softchk\|stamp\|element\|parasitic cap\|attribute cap\)"
"syn match   gmdstatement   "^ *\(flagnon45\|lextract\|equation\|lpeselect\|lpechk\|attach\)"
"syn match   gmdStatement   "\(temporary\|connect\)-layer"
"syn match   gmdStatement   "program-dir"
"syn match   gmdStatement   "status-command"
"syn match   gmdStatement   "batch-queue"
"syn match   gmdStatement   "cnames-csen"
"syn match   gmdStatement   "filter-lay-opt"
"syn match   gmdStatement   "filter-sch-opt"
"syn match   gmdStatement   "power-node"
"syn match   gmdStatement   "ground-node"
"syn match   gmdStatement   "subckt-name"

"syn match   gmdGsel		"\*description"
"syn match   gmdGsel		"\*input-layer"
"syn match   gmdGsel		"\*operation"
"syn match   gmdGsel		"\*end"

syn match   gmdComment "(\p*)"
syn match   gmdComment "^\s*# .*"


syn match   gmdPreProc "^\s*#\w.*"

"Modify the following as needed.  The trade-off is performance versus
"functionality.
syn sync lines=50

" Define the default highlighting.
" For version 5.7 and earlier: only when not done already
" For version 5.8 and later: only when an item doesn't have highlighting yet
if version >= 508 || !exists("did_gmd_syn_inits")
  if version < 508
    let did_gmd_syn_inits = 1
    command -nargs=+ HiLink hi link <args>
  else
    command -nargs=+ HiLink hi def link <args>
  endif

  HiLink gmdBlock      Special
  HiLink gmdIdentifier Identifier
  HiLink gmdComment    Comment
  HiLink gmdPreProc    PreProc
  HiLink gmdBoolean    Boolean
"  HiLink gmdAdd        Type
  HiLink gmdBCrea      String
  HiLink gmdSelected   Type
  HiLink gmdStatement  Statement
  HiLink gmdStop       Todo
  HiLink gmdGsel       Type

  delcommand HiLink
endif

let b:current_syntax = "gmd"
" vim: ts=8 
