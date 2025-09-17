# exportFASTA() on gene-sequence data.frame

    Code
      exportFASTA(test_df, stdout())
    Output
      >gene1
      ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
      
      >gene2
      GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
      

# exportFASTA() on 4 columns data

    Code
      exportFASTA(test_df, stdout())
    Output
      >chromosome_chr1_start_1_end_41
      ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
      
      >chromosome_chr1_start_55_end_95
      GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
      

