function Produce_list_HCP_subjects_codes_paper_GSR(Work_dir)

    fileID=fopen([Work_dir '/List_datasetHCP_subject_codes.txt'],'a');
    List_datasetHCP_subject_codes=   [100408 101107 101309 101410 101915 103515 103818 104416 104820 105115 105216 105620 ...
                             105923 107321 107422 107725 108525 109123 109325 110007 110613 111312 112112 112314 ...
                             112920 113215 113821 114419 115320 116221 116524 116726 117324 117930 118124 118225 ...
                             118528 118932 119126 119732 120111 120515 120717 121315 121416 121618 121921 122620 ...
                             122822 123824 124826 126628 128127 129028 129129 129533 129634 130013 130417 130619 ...
                             131217 131823 132017 133827 134021 134223 134324 134728 134829 136732 137633 138534 ...
                             138837 139233 140117 140319 141119 141826 142828 143426 144125 145127 145834 146432 ...
                             146533 146634 147030 150625 151425 151526 151627 151829 153025 153429 153833 154229 ...
                             154532 154734 155635 156233 156334 156435 156637 157942 158136 158843 159239 159340 ...
                             159845 160729 162733 163129 163432 164030 164131 164939 165234 165638 166640 167238 ...
                             167743 169343 169444 169949 170631 172534 172938 173435 173738 173940 175035 175338 ...
                             175742 176441 176542 176744 177241 177342 177746 178243 179346 180129 180937 183337 ...
                             185947 186444 189450 191841 192035 192136 192439 192540 192641 192843 193441 194140 ...
                             194645 194847 195041 196144 196750 197348 197651 198249 198451 198855 199453 199655 ...
                             200109 200210 200311 200917 201111 201717 202719 203923 204319 204420 204622 206222 ...
                             207123 208428 209935 210011 211215 211720 212015 212116 212217 212419 212823 213421 ...
                             214423 214524 214625 214726 220721 221218 233326 236130 250932 256540 268850 270332 ...
                             275645 280941 284646 285345 286650 289555 290136 295146 303119 303624 308331 310621 ...
                             316633 317332 321323 322224 334635 336841 341834 348545 361234 361941 377451 381038 ...
                             389357 395251 395756 397154 412528 414229 422632 424939 441939 449753 456346 459453 ...
                             465852 467351 475855 485757 512835 513736 517239 520228 521331 523032 529953 541943 ...
                             548250 552544 553344 555651 561444 562345 566454 567052 571144 573249 576255 579665 ...
                             579867 580044 581349 581450 583858 585862 588565 592455 597869 599469 599671 604537 ...
                             614439 616645 618952 620434 622236 626648 633847 638049 644044 645450 647858 656253 ...
                             656657 657659 660951 663755 664757 665254 671855 677968 679568 680957 695768 706040 ...
                             713239 724446 725751 727553 735148 748662 749058 753150 753251 756055 766563 767464 ...
                             769064 770352 771354 773257 782561 783462 788876 803240 810843 826353 833148 835657 ...
                             837560 841349 843151 844961 845458 849264 856463 861456 871762 871964 873968 889579 ...
                             894067 896778 901038 901139 904044 907656 910241 917255 919966 926862 942658 947668 ...
                             952863 955465 957974 958976 959574 966975 972566 984472 987983 991267 992673 993675 ...
                             996782];
    for subject=1:361
        fprintf(fileID,'%i \n',List_datasetHCP_subject_codes(subject));
    end
    fclose(fileID);
    
end