# Class 19
Harshitha Palacharla

``` r
library(bio3d)

sd <- read.fasta("A17775151_mutant_seq.fa")
sd
```

                   1        .         .         .         .         .         60 
    wt_healthy     MPPRPSSGELWGIHLMPPRILVECLLPNGMIVTLECLREATLITIKHELFKEARKYPLHQ
    mutant_tumor   MPPRPSSGELWGIHLMPPRILVECLLPNGMIVTLECLREATLITIKHELFKEARKYPLHQ
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     LLQDESSYIFVSVTQEAEREEFFDETRRLCDLRLFQPFLKVIEPVGNREEKILNREIGFA
    mutant_tumor   LLQDESSYIFVSVTQEAEREEFFDETRRLCDLRLFQPFLKVIEPVGNREEKILNREIGFA
                   ************************************************************ 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     IGMPVCEFDMVKDPEVQDFRRNILNVCKEAVDLRDLNSPHSRAMYVYPPNVESSPELPKH
    mutant_tumor   IGMPVCEFDMVKDPEVQDFRRNILNVCKEAVDLRDLNSPHSRAMYVYPPNVESSPELPKH
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     IYNKLDKGQIIVVIWVIVSPNNDKQKYTLKINHDCVPEQVIAEAIRKKTRSMLLSSEQLK
    mutant_tumor   IYNKLDKGQIIVVIWVIVSPNNDKQKYTLKINHDCVPEQVIAEAIRKKTRSMLLSSEQLK
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     LCVLEYQGKYILKVCGCDEYFLEKYPLSQYKYIRSCIMLGRMPNLMLMAKESLYSQLPMD
    mutant_tumor   LCVLEYQGKYILKVCGCDEYFLEKYPLSQYKYIRSCIMLGRMPNLMLMAKESLYSQLPMD
                   ************************************************************ 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     CFTMPSYSRRISTATPYMNGETSTKSLWVINSALRIKILCATYVNVNIRDIDKIYVRTGI
    mutant_tumor   CFTMPSYSRRISTATPYMNGETSTKSLWVINSALRIKILCATYVNVNIRDIDKIYVRTGI
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     YHGGEPLCDNVNTQRVPCSNPRWNEWLNYDIYIPDLPRAARLCLSICSVKGRKGAKEEHC
    mutant_tumor   YHGGEPLCDNVNTQRVPCSNPRWNEWLNYDIYIPDLPRAARLCLSICSVKGRKGAKEEHC
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .         .         .         480 
    wt_healthy     PLAWGNINLFDYTDTLVSGKMALNLWPVPHGLEDLLNPIGVTGSNPNKETPCLELEFDWF
    mutant_tumor   PLAWGNINLFDYTDTLVSGKMALNLWPVPHGLEDLLNPIGVTGSNPNKETPCLELEFDWF
                   ************************************************************ 
                 421        .         .         .         .         .         480 

                 481        .         .         .         .         .         540 
    wt_healthy     SSVVKFPDMSVIEEHANWSVSREAGFSYSHAGLSNRLARDNELRENDKEQLKAISTRDPL
    mutant_tumor   SSVVKFPDMSVIEEHANWSVSREAGFSYSHAGLSNRLARDNELRENDKEQLKAISTRDPL
                   ************************************************************ 
                 481        .         .         .         .         .         540 

                 541        .         .         .         .         .         600 
    wt_healthy     SEITEQEKDFLWSHRHYCVTIPEILPKLLLSVKWNSRDEVAQMYCLVKDWPPIKPEQAME
    mutant_tumor   SEITEQEKDFLWVHRHYCVTIPEILPKLLLSVEWNSRDEVAQMYCLVKDWPPIKPEQAME
                   ************ ******************* *************************** 
                 541        .         .         .         .         .         600 

                 601        .         .         .         .         .         660 
    wt_healthy     LLDCNYPDPMVRGFAVRCLEKYLTDDKLSQYLIQLVQVLKYEQYLDNLLVRFLLKKALTN
    mutant_tumor   LLDCNYPDPMVRGFARRCLEKYLTDDKLSQYLIQLVQVLKYYQYLDNLLVRFLLKKALTN
                   *************** ************************* ****************** 
                 601        .         .         .         .         .         660 

                 661        .         .         .         .         .         720 
    wt_healthy     QRIGHFFFWHLKSEMHNKTVSQRFGLLLESYCRACGMYLKHLNRQVEAMEKLINLTDILK
    mutant_tumor   QRIGHFFFWHLKSEMHNKTVSQRFGLLLESYCRACGMYLKHLNRQVEAMEKLINLTDILK
                   ************************************************************ 
                 661        .         .         .         .         .         720 

                 721        .         .         .         .         .         780 
    wt_healthy     QEKKDETQKVQMKFLVEQMRRPDFMDALQGFLSPLNPAHQLGNLRLEECRIMSSAKRPLW
    mutant_tumor   QEKKDETQKVQMKFLVEQMRRPDFMDALQGFLSPLNPAHQLGNLRLEECRIMSSAKRPLW
                   ************************************************************ 
                 721        .         .         .         .         .         780 

                 781        .         .         .         .         .         840 
    wt_healthy     LNWENPDIMSELLFQNNEIIFKNGDDLRQDMLTLQIIRIMENIWQNQGLDLRMLPYGCLS
    mutant_tumor   LNWENPDIMSELLFQNNEIIFKNGDDLRQDMLTLQIIRIMENIWQNQGLDLRMLPYGCLS
                   ************************************************************ 
                 781        .         .         .         .         .         840 

                 841        .         .         .         .         .         900 
    wt_healthy     IGDCVGLIEVVRNSHTIMQIQCKGGLKGALQFNSHTLHQWLKDKNKGEIYDAAIDLFTRS
    mutant_tumor   IGDCVGLIEVVRNSHTIMQIQCKGGLKGALQFNSHTLHQWLKDKNKGEIYDAAIDLFTRS
                   ************************************************************ 
                 841        .         .         .         .         .         900 

                 901        .         .         .         .         .         960 
    wt_healthy     CAGYCVATFILGIGDRHNSNIMVKDDGQLFHIDFGHFLDHKKKKFGYKRERVPFVLTQDF
    mutant_tumor   CAGYCVATFILGIGDRHNSNIMVKDDGQLFHIDFGHFLDHKKKKFGYKRERVPFVLTQDF
                   ************************************************************ 
                 901        .         .         .         .         .         960 

                 961        .         .         .         .         .         1020 
    wt_healthy     LIVISKGAQECTKTREFERFQEMCYKAYLAIRQHANLFINLFSMMLGSGMPELQSFDDIA
    mutant_tumor   LIVISKGAQECTKTREFERFQEMCYKAYLAIRQHANLFINLFSMMLGSGMPELQSFDDIA
                   ************************************************************ 
                 961        .         .         .         .         .         1020 

                1021        .         .         .         .       1068 
    wt_healthy     YIRKTLALDKTEQEALEYFMKQMNDAHHGGWTTKMDWIFHTIKQHALN
    mutant_tumor   YIRKTLALDKTEQEALEYFMKQMNDAHHGGWTTKMDWIFHTIKQHALN
                   ************************************************ 
                1021        .         .         .         .       1068 

    Call:
      read.fasta(file = "A17775151_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 1068 position columns (1068 non-gap, 0 gap) 

    + attr: id, ali, call
