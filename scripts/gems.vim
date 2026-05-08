" Vim syntax file
" Language:	GEMS
" Maintainer:	alexis paz <apaz@fcq.unc.edu.ar>
" Extensions:   *.gms
" Comment:      GEMS is a input script language for the General Molecular Simulator code.

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

" Ignore case
syn case ignore

" A bunch of useful gms keywords

"syn match   gmsIdentifier

syn keyword gmsBlock       bloque 
syn keyword gmsBlock       fin
syn keyword gmsBlock       repeat

syn keyword gmsStop        stop

syn keyword gmsIdentifier  out outfile box neigh time
syn keyword gmsSelected    set interact under evolve
syn match   gmsAdd         "\p* add"
syn match   gmsSelected    "\s<\s"

" Math and variables
" syn match   gmsBCrea     "{\p\{-}}"
syn match   gmsvar	  "\$[^ }]*"
syn match   gmsvar	  "\$[^ }]*"
syn region  gmsmath         start='{' end='}' contains=gmsvar

"syn keyword gmsGsel        sys group 
syn match   gmsBCrea       "+\s*\p*"
syn match   gmsBCrea       "-\s*\p*"
syn match   gmsGsel        "group\b*\d"
syn match   gmsGsel        "^\s*>\s*\p*" contains=gmsBCrea
syn match   gmsGsel        "^\s*\^\s*\p*" contains=gmsBCrea
syn match   gmsGsel        "^\s*>>\s*\p*"
syn match   gmsGsel        "^\s*\^\~\s*\p*"

" \zs and \ze set the start and the final of the submatch to higlight
syn match   gmsSelected    "\d\s\s*\zs<>\ze\s"
syn match   gmsSelected    "\d\s\s*\zs>\ze\s"

syn keyword gmsStatement   partemp dinamica hd lbfgs cg lp

"syn match gmsBoolean	"\ \s*\(T\|F\)\s* "
"
"
"syn match   gmsStatement   "flag-\(non45\|acuteangle\|offgrid\)"
"syn match   gmsStatement   "text-pri-only"
"syn match   gmsStatement   "[=&]"
"syn match   gmsStatement   "\[[^,]*\]"
"syn match   gmsstatement   "^ *\(sel\|width\|ext\|enc\|area\|shrink\|grow\|length\)"
"syn match   gmsstatement   "^ *\(or\|not\|and\|select\|size\|connect\|sconnect\|int\)"
"syn match   gmsstatement   "^ *\(softchk\|stamp\|element\|parasitic cap\|attribute cap\)"
"syn match   gmsstatement   "^ *\(flagnon45\|lextract\|equation\|lpeselect\|lpechk\|attach\)"
"syn match   gmsStatement   "\(temporary\|connect\)-layer"
"syn match   gmsStatement   "program-dir"
"syn match   gmsStatement   "status-command"
"syn match   gmsStatement   "batch-queue"
"syn match   gmsStatement   "cnames-csen"
"syn match   gmsStatement   "filter-lay-opt"
"syn match   gmsStatement   "filter-sch-opt"
"syn match   gmsStatement   "power-node"
"syn match   gmsStatement   "ground-node"
"syn match   gmsStatement   "subckt-name"

"syn match   gmsGsel		"\*description"
"syn match   gmsGsel		"\*input-layer"
"syn match   gmsGsel		"\*operation"
"syn match   gmsGsel		"\*end"

syn match   gmsComment "(\p*)"
syn match   gmsComment "^\s*# .*"


syn match   gmsPreProc "^\s*#\w.*"

"Modify the following as needed.  The trade-off is performance versus
"functionality.
syn sync lines=50

" Define the default highlighting.
" For version 5.7 and earlier: only when not done already
" For version 5.8 and later: only when an item doesn't have highlighting yet
if version >= 508 || !exists("did_gms_syn_inits")
  if version < 508
    let did_gms_syn_inits = 1
    command -nargs=+ HiLink hi link <args>
  else
    command -nargs=+ HiLink hi def link <args>
  endif

  HiLink gmsBlock      Special
  HiLink gmsIdentifier Identifier
  HiLink gmsComment    Comment
  HiLink gmsPreProc    PreProc
  HiLink gmsBoolean    Boolean
"  HiLink gmsAdd        Type
  HiLink gmsmath          String
  HiLink gmsBCrea      String
  HiLink gmsvar        Type
  HiLink gmsSelected   Type
  HiLink gmsStatement  Statement
  HiLink gmsStop       Todo
  HiLink gmsGsel       Type

  delcommand HiLink
endif

let b:current_syntax = "gms"
" vim: ts=8 
