>> How to use shaep <<

Takeshi Kawabata
2013年  5月 30日 木曜日

Usage examples

ShaEP lists available options when called with '-h' or '--help'.

    Superimpose A on B

    Here B is the reference (query) molecule that does not move. Partial atomic charges should be present in the file. Resulting structures will be written into 'out.sdf' (in SD format) and the similarity index to 'similarity.txt'

    shaep -q B.mol2 A.mol2 -s out.sdf similarity.txt

    Now the same, but using only volume overlap optimization, which does not require the partial atomic charges:

    shaep --onlyshape -q B.mol2 A.mol2 -s out.sdf shapesimilarity.txt


>> 典型的な実行コマンド　<<

shaep -q [query(reference)_molecular_file] [target_molecular_file] -s [superimposed target file in SDF] --output-file [out]


  *query  : template (fixed pose and position)
  *target : to be transformed

query, targetともに、PDB, MOL2, SDFを指定可能。target_molecular_fileとして、複数の配座が入ったmol2ファイルを
指定した場合、それらすべてが読み込み、計算される。

・-sオプションで出力されるsuperimposeされた分子構造は、必ずSDF形式で出力される。ファイルフォーマットの変更はできないらしい。

>>他の有用なオプション<<

 --nStructures arg              The number of superimpositions to output for 
                                 each molecule. Note that the structures are 
                                 output in the order of increasing similarity 
                                 index, i.e., worst first. Zero imposes no 
                                 limit. Default = 1

　出力される重ね合わせ構造の数。しかし、説明をよく読むと、出力はsimilarity indexの悪い順だと
　書いてある。それでは、意味がないので、-nStructures 0 として、全部のconformersを出力し、
　最後に出力されたconformerを最良構造として用いることにする。


>>注意点<<

・[query_molecule_file]として、PDBファイルを使う場合、いくつか注意が必要。

　(1) ATOM/HETATM行の末尾にある元素記号が必要。末尾の元素記号がないと、

　　Exception while parsing PDB file on line 1:
　　: Line is truncated.

　というエラーを出力して、不正終了してしまう。

  (2) ATOM/HETATM行に、'A','B',..などのalt_idがあった場合、例えば、2cchのATPのように、
　　　
HETATM 9347  PG AATP A1297      19.607   9.130   6.848  0.50 74.58           P
HETATM 9351  PB AATP A1297      18.847   6.517   6.137  0.50 74.02           P
HETATM 9355  PA AATP A1297      19.100   4.556   8.185  0.50 70.82           P
HETATM 9356  O1AAATP A1297      20.366   3.802   8.516  0.50 70.49           O
HETATM 9357  O2AAATP A1297      18.102   4.850   9.282  0.50 70.82           O
HETATM 9358  O3AAATP A1297      19.514   5.957   7.498  0.50 72.44           O
　　　を用いた場合、

terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check
アボートしました (コアダンプ)
　　　というエラーを出力して不正終了してしまう。

　
　(3) HEADER行が短すぎる場合、エラーが表示される。
　　例えば、
　　HEADER    COMPLEX (TRANSFERASE/CYCLIN)            14-JUL-96   1FIN
　　では正常終了するが、
　　HEADER    COMPLEX (TRANSFERASE/CYCLIN)
　　となっていると、
　　　　Exception while parsing PDB file on line 2:
　　　　: Line is truncated.
　　というエラーが表示される。ただし、エラーが表示されても計算は正常に終了する。


・オプション --output-file [out]を指定すると、[out]というファイルと[out]_hits.txtというファイルの二つの同一のファイルが
　なぜか出力される。



・-sオプションで出力されるSDFファイルには、以下のような、Similarityのタグが付加されている。

 11 38  1  0  0  0  0
 11 39  1  0  0  0  0
M  END
>  <Similarity_ESP>
0.37373834481693136

>  <Similarity_best>
0.52712155649975401

>  <Similarity_hit>
1NNC

>  <Similarity_shape>
0.68050476818257655

・複数のconformerが出力される場合、そのSimilarityは以下のようになる。つまり、
　(1) Similarity_bestは、Similarity_ESPとSimilarity_shapeの平均である。
　(2) Similarity_bestの低い順（悪い順）に出力される。

>  <Similarity_ESP> 0.37373834481693136 >  <Similarity_best> 0.52712155649975401 >  <Similarity_shape> 0.68050476818257655
>  <Similarity_ESP> 0.3633378003765621  >  <Similarity_best> 0.53321223751007796 >  <Similarity_shape> 0.70308667464359387
>  <Similarity_ESP> 0.29911840798627609 >  <Similarity_best> 0.53350089244962673 >  <Similarity_shape> 0.76788337691297726
>  <Similarity_ESP> 0.3349733990299763  >  <Similarity_best> 0.53778407205260292 >  <Similarity_shape> 0.74059474507522949
>  <Similarity_ESP> 0.35206549340203175 >  <Similarity_best> 0.55413884280789394 >  <Similarity_shape> 0.75621219221375624
:
>  <Similarity_ESP> 0.41196966301816823 >  <Similarity_best> 0.61964514274686389 >  <Similarity_shape> 0.8273206224755596



