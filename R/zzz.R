set_relation_dissimilarity_method("symdiff",
                                  "symmetric difference distance",
                                  .relation_dissimilarity_symdiff_method,
                                  .relation_dissimilarity_symdiff_xtrafo,
                                  .relation_dissimilarity_symdiff_xtrafo,
                                  .relation_dissimilarity_symdiff_params)

set_relation_dissimilarity_method("SD",
                                  "symmetric difference distance",
                                  .relation_dissimilarity_symdiff_method,
                                  .relation_dissimilarity_symdiff_xtrafo,
                                  .relation_dissimilarity_symdiff_xtrafo,
                                  .relation_dissimilarity_symdiff_params)

set_relation_dissimilarity_method("CKS",
                                  "Cook-Kress-Seiford distance",
                                  .relation_dissimilarity_CKS_method,
                                  .relation_dissimilarity_CKS_xtrafo,
                                  .relation_dissimilarity_CKS_xtrafo)

set_relation_dissimilarity_method("CS",
                                  "Cook-Seiford distance",
                                  .relation_dissimilarity_CS_method,         
                                  .relation_dissimilarity_CS_xtrafo,
                                  .relation_dissimilarity_CS_xtrafo)

set_relation_dissimilarity_method("score",
                                  "score-based distance",
                                  .relation_dissimilarity_score_method,
                                  .relation_dissimilarity_score_xtrafo,
                                  .relation_dissimilarity_score_xtrafo,
                                  .relation_dissimilarity_score_params)

set_relation_dissimilarity_method("manhattan",
                                  "Manhattan distance",
                                  .relation_dissimilarity_manhattan_method,
                                  relation_incidence,
                                  relation_incidence,
                                  .relation_dissimilarity_manhattan_params)

set_relation_dissimilarity_method("euclidean",
                                  "Euclidean distance",
                                  .relation_dissimilarity_euclidean_method,
                                  relation_incidence,
                                  relation_incidence,
                                  .relation_dissimilarity_euclidean_params)

set_relation_dissimilarity_method("Jaccard",
                                  "Jaccard distance",
                                  .relation_dissimilarity_Jaccard_method,
                                  relation_incidence,
                                  relation_incidence,
                                  .relation_dissimilarity_Jaccard_params)

set_relation_dissimilarity_method("PC",
                                  "paired comparison distance",
                                  .relation_dissimilarity_PC_method,
                                  .relation_dissimilarity_PC_xtrafo,
                                  .relation_dissimilarity_PC_ytrafo,
                                  .relation_dissimilarity_PC_params)
