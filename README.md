# seqFF
计算NIPT中的胎儿分数（fetal fraction）
脚本基于Sung K. Kim et al.的
[Determination of fetal DNA fraction from the plasma of pregnant women using sequence read counts](https://doi.org/10.1002/pd.4615)


## 依赖以下R包和软件
1. tidyverse
2. stringr
3. GenomicRanges
4. tools
5. magrittr
6. samtools

## 输入
一个bam文件，建议去重

## 输出
一个文本文件，一行四列，第一列为bam文件名，第二，三列分别为Elastic Net和WRSC拟合后的胎儿分数，第四列为两者的平均值，可用作最终结果

## 使用方法
单端测序
./seqFF.R SE [bam_file]

双端测序
./seqFF.R PE [bam_file]
