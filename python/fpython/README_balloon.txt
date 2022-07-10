>> How to use balloon  <<
2013年  5月 13日 月曜日
Takeshi Kawabata,

** 典型的なコマンド例 **

balloon --input-file [*.sdf|*.mol2|*.smi|*.pdb] --output-file [*.sdf|*.mol2|*.smi]


** Usefull options **

 
 -f [ --forcefield ] arg     A filename to read for MMFF94 force field 
                              parameters. Alternatively, you can set 
                              environment variable named BALLOON_FORCEFIELD to 
                              point to the parameters file. The command line 
                              option overrides the environment variable.
　MMFF94のエネルギーパラメータファイルの所在。"-f ~/tool/balloon/MMFF94.mff"などと設定する。
　あるいは、環境変数BALLOON_FORCEFIELDに設定してもいいらしい。

 --nGenerations arg             GA: The maximum number of generations. Default= 20
　　　　　　　　　　　　　　　　　　　GAの世代数。デフォルトは20だが、ShaEPの論文では200に設定したとある。

 -c [ --nconfs ] arg (=1)       Number of conformers to generate or the 
                                 initial population size if using GA. Zero will
                                 cause the input structure to be written out 
                                 with partial atomic charges and energy as 
                                 calculated by the MMFF94-like force field.
　　　　　　　　　　　　　　　　　出力するconformerの数。ShaEPの論文では20に設定したとある。
                               ただし、実際に計算してみると、指定した-cの値と少し異なる数の構造が生成される場合が結構ある...
                               例えば、PDBのN20の場合、-c 10でも８個、-c 20でも25個数生成される。

 -i [ --maxiter ] arg           Maximum allowed number of iterations for 
                                 conjugate gradient structure optimization in 
                                 the template geometry generation. Default = 100
　　　　　　　　　　　　　　　　　　共役勾配法の最大反復数。デフォルトは100だが、ShaEPの論文では20に設定したとある。

　--singleconf                   Output only the lowest-energy conformation, 
                                regardless of the population size.

 -H [ --adjusthydrogens ]       Add hydrogens to the structures according to 
                                 the octet rule. Hydrogens are always added to 
                                 structures parsed from SMILES.

   --randomSeed arg               Seed the pseudo-random number generator. Range
                                 [1, 4294967295). Default value taken from  clock.
　　　　　　　　　　　　　　　　　　　乱数のシード。デフォルトは時計から取得。再現性を維持するためには、
　　　　　　　　　　　　　　　　　　　適当な数で固定したほうがいいかもしれない。

キラリティをフリップさせる機能もついている。が、デフォルトの確率は0.05になっており、
それほど高くはない。多くのエナンチオマーを含む化合物の場合、この確率を上げたほうがいいかも
しれない。

  --pStereoMutation arg          GA: Mutation probability for inverting a 
                                 stereochemical center. Default =   0.050000000000000003

  --pPyramidMutation arg         GA: Mutation probability for inverting a 
                                 pyramidal center. Default = 0.050000000000000003

  --pRingFlipMutation arg        GA: Mutation probability for a ring flip. 
                                 Default = 0.050000000000000003


