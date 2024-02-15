This file explain what's included in the sheets of "gatk_output.xlsx" file.

1-	All_variants: There are 183 variants and 569 RI(AI)Ls.

2-	rm_sites_Ls: 14 variants are excluded from the rest of analysis (Table 1).

Table 1. List of the variants that are removed from the rest of analysis:
#CHROM	POS
I	957607
I	14170805
II	4078227
III	5985312
III	10819732
IV	3359756
IV	3378546
IV	11312467
V	9033433
V	13511859
X	6064683
X	9735171
X	12722852
X	14067695

3-	rm_RI(AI)Ls: The RILs with the number of missing genotypes ≥ 160 were removed from the future analysis,23 RILs (Table 2).

Table 2. The lines with ≥ 160 missing variants.
No	Lines
1	RIL_104_S239
2	RIL_111_S73
3	RIL_112_S109
4	RIL_186_S31
5	RIL_195_S61
6	RIL_203_S244
7	RIL_21_S312
8	RIL_270_S89
9	RIL_294_S228
10	RIL_301_S250
11	RIL_316_S341
12	RIL_326_S330
13	RIL_335_S333
14	RIL_341_S323
15	RIL_347_S191
16	RIL_35_S38
17	RIL_362_S155
18	RIL_39_S200
19	RIL_402_S343
20	RIL_407_S365
21	RIL_41_S68
22	RIL_42_S103
23	RIL_62_S6

4-	GT_extraction: The genotypes were extracted here. There are 546 RI(AI)Ls and 169 variants.

5-	Final: The final data including all the variants. There are 27 lines that have GT information but there are no fitness data for them (Table 3). In this sheet the wild type = 0, mutatnt = 1 and missing data and het sites = NA

Table 3. For these lines there are GT information and there are not any fitness data.
No	Line
1	RIL_219
2	RIL_411
3	RIL_412
4	RIL_413
5	RIL_416
6	RIL_417
7	RIL_418
8	RIL_419
9	RIL_421
10	RIL_422
11	RIL_423
12	RIL_424
13	RIL_426
14	RIL_427
15	RIL_428
16	RIL_429
17	RIL_430
18	RIL_431
19	RIL_432
20	RIL_433
21	RIL_434
22	RIL_435
23	RIL_436
24	RIL_437
25	RIL_438
26	RIL_439
27	RIL_440

6-	SNPs: This sheet includes the base-substitutions only.

7-	INDELs: This sheet includes the indels only.

8-	Missing_No: This sheet has the information about missing data.
