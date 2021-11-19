# Vseg Pipeline

## Data Download
Test data could be downloaded at https://drive.google.com/drive/folders/1Ecn0xosDqCvQwRZSXn7aHEIOmWQm9wC6?usp=sharing. there are two files:

### OB36suniq.coord
it is a barcode list, each row represents one barcode, each row include barcodes' info as below:

| index | CoordX  | CoordY  | UMIcounts  | mtCounts  |  intronCounts  |
|---|---|---|---|---|---|
| 873491 | 24000.0115 | 6118.725 | 4 | 0 | 1 |
| 873492 | 24000.028 | 4329.86 | 1 | 1 | 0 |
| 873493 | 24000.0385 | 3805.31 | 2 | 0 | 0 |
| 873494 | 24000.05 | 9380.65 | 1 | 1 | 0 |

### OB36_sub.final
it is a transcript reads list, each row represent one transcript, it includes transcripts' info as below:

|Fov row|Fov column|X coordinate|Y coordinate|spatical barcode|barcode|Forward(0) or reverse(16) match|No. Chrosome|gene start location|gene matching pattern| |gene id|gene name|gene type|intro ratio|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|014|004|29126.2|6833.4|TATTTCC|AATTAGCCACGAACATGAGGTCAT|0|1|39370082|58M814N92M|1|ENSMUSG00000073702|Rpl31|protein_coding|0.37|
|014|004|29515.71|6708|ACATTGCGG|GTTACTGCATAAATTTGTGTCTGC|16|9|109078311|115M3576N35M|1|ENSMUSG00000091537|Tma7|protein_coding|1.00|
|013|004|29624.67|6888.33|TTGATGC|ATGAATTGCTTTGTTACGTAAAGT|16|12|111961438|101M894N24M773N25M|1|ENSMUSG00000021290|Atp5mpl|protein_coding|0.00|
  
