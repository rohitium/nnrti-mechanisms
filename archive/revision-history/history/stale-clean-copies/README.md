# Stale CLEAN copies — do not edit, do not submit

Moved aside 26-Oct-2026. Nobody is working on them and keeping two copies in
step was costing more than it was worth.

They are **out of date** and already differed from the marked copies in nine
places before this revision began (British/American spelling in the
convergence-diagnostics paragraph, and different wording in Statistical
analysis and in "Several genotypes show inconsistencies"). The SI copies also
diverge in convention: the CLEAN SI numbers its section headings, the marked
one does not.

**Regenerate them from the marked copies at the end**, by stripping the red
character colour, rather than editing them in parallel. That is one pass over
`word/document.xml` removing `<w:color w:val="FF0000"/>`, and it guarantees the
two match word for word — which is what the submission requires and what
maintaining them by hand kept failing to deliver.
