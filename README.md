# TODO

- incorporate txdb tbls in src_organism

- ?? collapse txchrom, txstart, txend, txstrand to tx column with
  format txchrom:txstart-txend:txstrand. Likewise for ex, cds

# DONE

-  all org views available in a single sql 'database' -- ATTACH one db
   to another?

- Update SQL tables

- Code cleanup

   - .get_tbl -- character(1); no need to check if view exists
   - documentation -- @rdname, etc.
   - tbl_org class also applies to tbl_txdb so...
