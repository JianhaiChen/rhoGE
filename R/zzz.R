## data.table's non-standard evaluation looks like undefined globals to R CMD check
utils::globalVariables(c(".", "BETA", "SE", "P", "SNP", "EA", "OA", "N", "gene", "tissue", "n_iv",
                         "b1", "b2", "se1", "se2", "p1", "p2", "pmin_", "p1_acat", "p2_acat",
                         "q1", "q2", "class", "type", "share", "z1", "z2", "chr", "pos", "block"))
