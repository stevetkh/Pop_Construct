  #!/bin/bash        

    for i in {2..32}
    do
      echo "Rscript BF_both_sexes_thiele_oag_single_year_spline_RW_bash_1_and_5.r $i" #check call to r script is as expected
      Rscript C:\\Users\\ktang3\\Desktop\\Imperial\\Pop_Construct\\BF_both_sexes_thiele_oag_single_year_spline_RW_bash_1_and_5.r $i
    done