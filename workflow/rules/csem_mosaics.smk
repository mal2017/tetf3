rule csem_peaks_regioner:
    """
    currently pulls directly from amarel, uses DSI library for pan
    but can use aggregated peaks from all libraries too.
    gro negcon is currently derived from all gro libraries

    also pulls "~/amarel-matt/tetf/subworkflows/tetf_tidal/results/beds/dgrp_tidal_insertions.unique.bb"
    from amarel
    """
    output:
        rds = "results/csem_mosaics/regioner.rds"
    script:
        "../scripts/csem_mosaics/regioner.R"