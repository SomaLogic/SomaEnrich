# SomaEnrich 0.1.0

* Initial package creation.
* `somaORA()` and `somaPrGSEA()` now return S3 objects of class `somaORA`
  and `somaPrGSEA`, respectively.
* `plotBubble()` is an S3 generic: ORA and GSEA results both use a ranked
  pathway bubble plot (fold enrichment or NES on the x-axis; bubble size is
  `% Overlap` or leading-edge %, respectively). GSEA plots are split into
  negative and positive NES panels, each sorted by adjusted p-value. The
  first argument was renamed from `ora_results` to `x`.
