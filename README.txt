###################################################################################################################################################################

Title: Hypoxia-induced changes in intragenic DNA methylation correlates with alternative splicing in breast cancer

List of commands used for various analysis in this paper.

Created:May the 4th of 2019
###################################################################################################################################################################
########################################################1.Differential Gene Expression Analysis######################################################################
###################################################################################################################################################################
mkdir Project
cd Project

mkdir 1.Diff_gene_exp_analysis
cd 1.Diff_gene_exp_analysis
#For getting most representative probe sets for Affymetrix U133 Plus 2 using Jetset package

#The file all_DRG.txt has been generated after TAC analysis and contains all differentially regulated genes in Microarray dataset GEO41491

#Removing duplicate rows from all_DRG.txt so that each diff. regulated gene appear once in the list
R
a=read.table(file="all_DRG.txt",header=T,sep="\t",quote="")
library(dplyr)
b=a%>%distinct(Gene.Symbol,.keep_all=T)
write.table(b,file="all_DRG_unique",sep="\t",append=F,row.names=F,quote=F)
q()
n
grep -v "Normoxia" all_DRG_unique |cut -f 7>all_DRG_gene_symbols
R
a=read.table(file="all_DRG_gene_symbols",header=F,as.is=TRUE)
library(jetset)
b=jmap('hgu133plus2',symbol=a$V1)
write.table(b,file="probeID_best",sep="\t",append=F,row.names=F,quote=F)
q()
n
R
a=read.table(file="probeID_best",header=T,sep="\t")
b=a%>%distinct(x,.keep_all=T)
write.table(b,file="probeID_best_uniq",sep="\t",row.names=F,quote=F)
q()
n
R
a=read.table(file="all_DRG.txt",sep="\t",header=T,quote="")
b=read.table(file="probeID_best_uniq",header=F,sep="\t")
c=merge(a,b,by.x="ID",by.y="V1")
d=c%>%distinct(ID,.keep_all=T)
write.table(d,file="all_DRG_m_probe_best_unique",sep="\t",row.names=F,quote=F)
q()
n
sed -i 's/ /_/g' all_DRG_m_probe_best_unique

#For getting the list of upregulated genes out of all differentially regulated genes
awk '$4<0{print $0}' all_DRG_m_probe_best_unique|sed 's/ /\t/g'> Upregulated_genes

cd ..

###################################################################################################################################################################
########################################################2.Identification of HIF-1α binding sites using ChIP-Seq data################################################
###################################################################################################################################################################


mkdir 2.Chip-seq_peak_calling
cd 2.Chip-seq_peak_calling
#For getting microarray probe coordinates (GRCh38) 
R
library('limma')
library(affy)
require(GEOquery)
library(Biobase)
library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
read.table(file="affy_list.txt",header=FALSE,as.is=TRUE)->J
x=getBM(attributes = c('affy_hg_u133_plus_2','hgnc_symbol','chromosome_name','start_position','end_position','band'),filters = 'affy_hg_u133_plus_2',values = (J$V1),mart=ensembl)
write.table(file="affy_hg_u133_plus_2.bed",x,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
q()
n

cp /home/erpl/Project/1.Diff_gene_exp_analysis/all_DRG_unique ../2.Chip-seq_peak_calling/

#merging probe coordinate file with all DRG
cut -f 1 all_DRG_unique > all_DRG_ID
R
a=read.table(file="affy_hg_u133_plus_2.bed",header=T,sep="\t",quote="")
b=read.table(file="all_DRG_ID",header=T,sep="\t")
c=merge(a,b,by.x="affy_hg_u133_plus_2",by.y="ID")
write.table(c,file="probe_m_DRG",append=F,sep="\t",row.names=F,quote=F)
q()
n
#removing scaffolds, rows without gene names
grep -v "CHR" probe_m_DRG > probe_m_DRG_final
grep -v "affy" probe_m_DRG_final|awk '{print $3,$4,$5,$1,$2}'|sed 's/ /\t/g'|sort -k1,1 -k2n,2>probe_m_DRG_final.bed
grep -v "p" probe_m_DRG_final.bed|grep -v "q"|grep -v "GL000194">probe_m_DRG_final_1.bed


#Bowtie buidling index
cd /home/erpl/Bowtie_index 
bowtie -build Homo_sapiens.GRCh38.dna_sm.toplevel.fa bowtie23.01.19

#GSM700944_ChIPSeq_HIF1_alpha_sorted.txt and GSM700946_ChIPSeq_Pre_immune_Control_For_HIF1_HIF2_alpha_sorted.txt were downloaded from GEO id# GSE28352

#Conversion of ChIP Seq files to FASTQ format
cd /home/erpl/Project/2.Chip-seq_peak_calling/
cut -f 9,10 GSM700944_ChIPSeq_HIF1_alpha_sorted.txt|awk '{print "@"NR"\n"$1"\n""+\n"$2}' > GSM700944_ChIPSeq_HIF1_alpha.fq
cut -f 9,10 GSM700946_ChIPSeq_Pre_immune_Control_For_HIF1_HIF2_alpha_sorted.txt|awk '{print "@"NR"\n"$1"\n""+\n"$2}' > GSM700944_ChIPSeq_Pre_immune_Control_For_HIF1_HIF2_alpha.fq

#Mapping with bowtie
bowtie -S /home/erpl/Bowtie_index/bowtie23.01.19 GSM700944_ChIPSeq_HIF1_alpha.fq > mapping_bowtie_6.4.19.sam
samtools sort -n mapping_bowtie_6.4.19.sam -o mapping_bowtie_sorted_6.4.19.sam
bowtie -S /home/erpl/Bowtie_index/bowtie23.01.19 GSM700944_ChIPSeq_Pre_immune_Control_For_HIF1_HIF2_alpha.fq > mapping_bowtie_input_6.4.19.sam
samtools sort -n mapping_bowtie_input_6.4.19.sam -o mapping_bowtie_input_sorted_6.4.19.sam

#Peak Calling with MACS2
mkdir macs2
cd /home/erpl/Project/2.Chip-seq_peak_calling/macs2

macs2 callpeak –t /home/erpl/Project/2.Chip-seq_peak_calling/mapping_bowtie_6.4.19.sam /home/erpl/Project/2.Chip-seq_peak_calling/mapping_bowtie_input_6.4.19.sam -f SAM -n 8.4.19_HIF1A_peaks --outdir /home/erpl/Project/2.Chip-seq_peak_calling/macs2
grep -v "CHR" 8.4.19_HIF1A_peaks.narrowPeak>8.4.19_HIF1A_peaks_latest.narrowPeak

#assigning nearest macs2 peak to each of the genes
cp /home/erpl/Project/2.Chip-seq_peak_calling/probe_m_DRG_final_1.bed ../macs2/
bedtools closest -a probe_m_DRG_final_1.bed -b 8.4.19_HIF1A_peaks_latest.narrowPeak -d>closest_bed_macs2
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$16}' closest_bed_macs2|sed 's/ /\t/g'|sed '1i chr\tprobe_start\tprobe_end\tProbe_ID\tGene_name\tchr\tPeak_start\tPeak_end\tDist_from_probe_to_peak'>closest_bed_macs2_with_headers

#TSS_all_genes.txt contains the transcription start sites (TSS) for all the differentially expressed genes downloaded from BioMart on Ensembl genome browser 96 (GRCh38.p12)
awk '{print $3,$4,$5,$6}' TSS_all_genes.txt|sed 's/ /\t/g'>TSS_all_genes_latest.txt
R
a=read.table(file="closest_bed_macs2_with_headers",header=TRUE,sep="\t",quote="")
b=read.table(file="TSS_all_genes_latest.txt",header=T,sep="\t")
c=merge(a,b,by.x="Probe_ID",by.y="AFFY_HG_U133_Plus_2_probe")
write.table(c,file="closest_bed_macs2_m_biomart_TSS",sep="\t",append=F,row.names=F,quote=F)
q()
n
awk '{print $12,$13,$1,$8,$9,$11,$9-$11}' closest_bed_macs2_m_biomart_TSS|sed 's/ /\t/g'>distance_peaks_TSS_macs2
#Finding genes that are located upstream or downstream 2 kb of the TSS
awk '$8<=2000{print $0}' distance_peaks_TSS_macs2|awk '$8>0{print $0}'>final_genes_with_downstream_2kb_macs2
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_downstream_2kb_macs2
awk '$8>=-2000{print $0}' distance_peaks_TSS_macs2|awk '$8<0{print $0}'>final_genes_with_upstream_2kb_macs2
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_upstream_2kb_macs2
cat final_genes_with_downstream_2kb_macs2 final_genes_with_upstream_2kb_macs2|grep -v "Gene_name"|cut -f 1|sort|uniq>gene_names_2kb_macs2_uniq

cd ..
#Peak Calling with GEM

mkdir GEM
cd GEM

#command used to perform peak calling with GEM
#java -Xmx10G -jar gem.jar --d Read_Distribution_default.txt --expt mapping_bowtie_6.4.19.sam --ctrl mapping_bowtie_input_6.4.19.sam --f BED –outNP --out /home/erpl/Project/2.Chip-seq_peak_calling/GEM
cd /home/erpl/Project/2.Chip-seq_peak_calling/GEM

grep -v "CHR" GEM_HIF1A.GPS_events.narrowPeak|grep -v "chrK">GEM_HIF1A.GPS_events_latest.narrowPeak

#assigning nearest GEM peak to each of the genes
bedtools closest -a probe_m_DRG_final_1.bed -b GEM_HIF1A.GPS_events_latest.narrowPeak -d>closest_bed_GEM
awk '{print $1,$2,$3,$4,$5,$7,$8,$16}' closest_bed_GEM|sed 's/ /\t/g'|sed '1i chr\tprobe_start\tprobe_end\tProbe_ID\tGene_name\tPeak_start\tPeak_end\tDist_from_probe_to_peak'>closest_bed_GEM_with_headers
R
a=read.table(file="closest_bed_GEM_with_headers",header=TRUE,sep="\t",quote="")
b=read.table(file="TSS_all_genes_latest.txt",header=T,sep="\t")
c=merge(a,b,by.x="Probe_ID",by.y="AFFY_HG_U133_Plus_2_probe")
write.table(c,file="closest_bed_GEM_m_biomart_TSS",sep="\t",append=F,row.names=F,quote=F)
q()
n
awk '{print $10,$11,$1,$2,$6,$7,$9,$7-$9}' closest_bed_GEM_m_biomart_TSS|sed 's/ /\t/g'>distance_peaks_TSS_GEM
#Finding genes that are located upstream or downstream 2 kb of the TSS
awk '$8<=2000{print $0}' distance_peaks_TSS_GEM|awk '$8>0{print $0}'>final_genes_with_downstream_2kb_GEM
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_downstream_2kb_GEM
awk '$8>=-2000{print $0}' distance_peaks_TSS_GEM|awk '$8<0{print $0}'>final_genes_with_upstream_2kb_GEM
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_upstream_2kb_GEM
cat final_genes_with_downstream_2kb_GEM final_genes_with_upstream_2kb_GEM|grep -v "Gene_name"|cut -f 1|sort|uniq>gene_names_2kb_GEM_uniq

#File containing both upstream and downstream genes is final_genes_2kb_GEM stored in /home/erpl/Project/2.Chip-seq_peak_calling/GEM

mkdir CCAT
cd CCAT
#Peak Calling with CCAT
peakranger ccat –d mapping_bowtie_6.4.19.sam -c mapping_bowtie_input_6.4.19.sam --format SAM –o HIF1A_CCAT_12.4.19
sort -k1,1 -k2,2n HIF1A_CCAT_12.4.19.bed>HIF1A_CCAT_12.4.19_sorted.bed

#assigning nearest CCAT peak to each of the genes
bedtools closest -a probe_m_DRG_final_1.bed -b HIF1A_CCAT_12.4.19_sorted.bed -d>closest_bed_CCAT
awk '{print $1,$2,$3,$4,$5,$7,$8,$12}' closest_bed_CCAT|sed 's/ /\t/g'|sed '1i chr\tprobe_start\tprobe_end\tProbe_ID\tGene_name\tPeak_start\tPeak_end\tDist_from_probe_to_peak'>closest_bed_CCAT_with_headers
R
a=read.table(file="closest_bed_CCAT_with_headers",header=TRUE,sep="\t",quote="")
b=read.table(file="TSS_all_genes_latest.txt",header=T,sep="\t")
c=merge(a,b,by.x="Probe_ID",by.y="AFFY_HG_U133_Plus_2_probe")
write.table(c,file="closest_bed_CCAT_m_biomart_TSS",sep="\t",append=F,row.names=F,quote=F)
q()
n
awk '{print $10,$11,$1,$2,$6,$7,$9,$7-$9}' closest_bed_CCAT_m_biomart_TSS|sed 's/ /\t/g'>distance_peaks_TSS_CCAT

#Finding genes that are located upstream or downstream 2 kb of the TSS
awk '$8<=2000{print $0}' distance_peaks_TSS_CCAT|awk '$8>0{print $0}'>final_genes_with_downstream_2kb_CCAT
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_downstream_2kb_CCAT
awk '$8>=-2000{print $0}' distance_peaks_TSS_CCAT|awk '$8<0{print $0}'>final_genes_with_upstream_2kb_CCAT
sed -i '1i Gene_name\tTranscript_name\tProbe_ID\tchr\tPeak_start\tPeak_end\tTSS\tDist_from_TSS_to_peak' final_genes_with_upstream_2kb_CCAT
cat final_genes_with_downstream_2kb_CCAT final_genes_with_upstream_2kb_CCAT|grep -v "Gene_name"|cut -f 1|sort|uniq>gene_names_2kb_CCAT_uniq

cd ..

mkdir Fig2a
cd Fig2a
cp /home/erpl/Project/2.Chip-seq_peak_calling/macs2/gene_names_2kb_macs2_uniq ../Fig2a
cp /home/erpl/Project/2.Chip-seq_peak_calling/GEM/gene_names_2kb_GEM_uniq ../Fig2a
cp /home/erpl/Project/2.Chip-seq_peak_calling/CCAT/gene_names_2kb_CCAT_uniq ../Fig2a
cp /home/erpl/Project/1.Diff_gene_exp_analysis/Upregulated_genes ../Fig2a
#"gene_names_2kb_CCAT_uniq", "gene_names_2kb_macs2_uniq", "gene_names_2kb_GEM_uniq" and "Upregulated_genes" were used to make venn diagram on interactivenn.com, the genes which overlapped between these 4 files were downloaded and saved as custom_sig_before_cleanup3Jun.txt

###################################################################################################################################################################
#########################################################3. Generation of Custom Hypoxia Signature###################################################################
###################################################################################################################################################################
#heatmap 
#For merging of custom list to expression at each timepoint

mkdir 3.Generation_of_custom_hyp_signature
cd 3.Generation_of_custom_hyp_signature
cp /home/erpl/Project/2.Chip-seq_peak_calling/Fig2a/custom_sig_before_cleanup3Jun.txt ../3.Generation_of_custom_hyp_signature
R
a=read.table(file="exp_each_timepoint",header=T,sep="\t")
b=read.table(file="custom_sig_before_cleanup3Jun.txt",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom_m_exp_each_timepoint_3Jun",sep="\t",append=F,row.names=F,quote=F)
q()
n
#merging msigdb list with expression level at each timepoint
#msigdb signature was downloaded from http://software.broadinstitute.org/gsea/msigdb/cards/HALLMARK_HYPOXIA.html
sed 's/AK3L1/AK4/' msigdb_signature.txt>msigdb_signature_alt_name.txt
R
a=read.table(file="exp_each_timepoint",header=T,sep="\t")
b=read.table(file="msigdb_signature_alt_name.txt",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="exp_each_timepoint_m_msigdb",sep="\t",append=F,row.names=F,quote=F)
q()
n

#Heatmap was made with online tool Morpheus

###################################################################################################################################################################
######################################################4.Collection of TCGA data ##############################################################################
###################################################################################################################################################################

mkdir 4.Collection_of_TCGA_data
cd 4.Collection_of_TCGA_data
wget https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/HiSeqV2_exon.gz
wget https://toil.xenahubs.net/download/tcga_Kallisto_est_counts.gz
wget https://toil.xenahubs.net/download/tcga_RSEM_gene_fpkm.gz
wget https://pancanatlas.xenahubs.net/download/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz
wget https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz
cd ..

###################################################################################################################################################################
######################################################5. Classification of TCGA Samples##############################################################################
###################################################################################################################################################################

#Classification of tumor samples

wget https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix.gz

#Extracting sample IDs corresponding to primary tumor or metastatic
grep -v "Solid Tissue Normal" BRCA_clinicalMatrix|cut -f 1 >tumor_sample_IDs

#Extracting gene expression values of only breast cancer samples
head -1 tcga_RSEM_gene_fpkm> datafile
for i in `cat tumor_sample_IDs|grep -v "sample_ID"`
do
columnname=`echo $i`
est=`cat datafile | tr "\t" "\n"|awk -v tar=$columnname '$1==tar{print NR}'`
echo $est
done
|tr "\n" ","
cut -f 3933,6677,10377,6530,6674,9880,9038,10503,8970,8109,8493,10132,5166,9716,8914,1341,2128,8839,6313,2605,9691,6082,6150,8760,2034,5209,1034,8048,5037,5998,10324,1538,5572,5766,1973,2774,8391,7957,422,5163,8320,641,2231,3686,8808,9616,2022,9746,2957,9538,7847,1955,6711,1917,5000,4665,3871,947,7654,2718,6427,872,7574,10436,6608,2976,9770,4839,8596,1054,1778,2285,6056,9247,8187,5091,8872,9076,2206,5897,2782,9530,9567,2580,5582,9360,3788,10524,7517,6744,5930,8731,4959,7761,4034,8571,563,841,384,7615,2991,8898,7119,9585,5736,2802,10294,2748,5539,5906,9705,2181,8903,28,9753,6821,8147,2675,5512,1772,7862,4105,5589,2952,8245,7227,5151,379,1493,2208,7109,4221,419,7752,8780,5624,1842,4989,4431,3740,8199,6087,9955,4302,3436,8517,8555,7619,3870,1941,5219,6981,2479,6920,235,7768,1340,8890,6987,8610,4856,7519,6937,3116,614,7053,9042,7998,10015,3259,6125,9945,3193,1979,1408,1216,8102,8430,5477,7001,3186,6791,6835,10396,6604,1049,4834,7831,5809,2085,8809,9682,302,9548,7357,10373,3618,3418,7196,7003,882,3923,7704,5837,5102,2118,6641,2773,1940,1759,6685,5885,9748,6713,1687,2834,10501,6649,3691,7425,10346,9513,250,10304,703,1728,797,2513,2926,1366,2584,3543,10528,6293,3347,3766,738,2494,6166,4563,4458,6242,2808,8706,9398,2628,9589,1004,2178,3868,2546,4540,7448,475,3049,6463,7287,6512,9949,3767,2712,150,4160,10027,8672,9023,6083,10344,3713,4491,9468,9627,2843,1585,3651,1829,4813,3138,9543,2526,684,6534,2060,5861,1647,10024,9193,245,4500,1522,1521,2772,6470,5711,1779,2595,5577,3980,5879,2089,1429,2244,8957,2195,5913,6108,8311,9138,4195,5234,10002,6117,1224,9136,9082,2339,9339,6446,5669,1580,9564,8214,1518,1644,5704,2779,3998,5874,2127,2937,3289,3487,2029,10505,6935,6902,185,5700,1961,4435,7363,5462,1823,5495,3942,2356,6500,2968,10519,6710,1790,5561,8525,9286,7859,1826,6218,8513,4755,2284,4822,6195,2220,1107,7810,2920,2130,1784,8313,1664,6740,3710,10469,4158,5370,3744,5379,10402,6545,6132,2626,2417,9151,105,8858,2633,8929,9856,3839,4634,5029,3303,5596,6591,10443,3032,1970,9396,2687,2864,6391,4821,8593,5554,7200,8453,10215,7673,5257,9890,3958,8657,10374,1417,5643,7444,481,4055,600,1847,639,8268,4233,8027,9166,8031,691,6409,3468,7547,4619,3841,446,1289,1558,5372,8314,9153,9969,5435,8275,7499,1621,8352,3450,4463,634,8180,1595,7830,8099,10386,6667,6373,1688,2542,4662,7649,7951,5013,801,7524,7069,3391,6015,6179,7169,468,10500,2831,3675,7890,1120,8478,4938,2364,4566,10280,8279,10314,1966,9611,404,299,7027,73,9840,6841,3273,8925,8076,689,3161,10098,4608,6974,7319,394,4204,9974,3503,545,3215,9673,5756,9591,350,8723,2070,9918,2478,4259,10048,247,5122,3887,5221,1461,8451,4694,2141,9745,1084,654,2688,8353,3743,10261,5484,9259,1093,10413,2699,4050,7962,2856,7093,5012,1203,559,6891,180,9930,8636,5105,4913,1108,3940,2058,503,169,3540,2094,9187,587,1883,7807,6663,9649,10406,1189,4104,1155,9002,2939,9577,3989,6021,10089,3279,7302,6704,557,493,4307,8491,2286,5276,8073,6068,8770,8973,1509,4680,489,8036,4495,7077,5389,3163,2436,814,1601,4478,8197,4446,7168,3369,9965,10225,10297,6421,6016,5093,2637,6404,10173,3335,6343,10498,2435,6216,3774,215,2514,6285,9270,10451,3706,675,7439,7081,3537,3493,3494,5265,5963,6266,7359,6067,2553,2350,1968,1476,5281,4582,8322,5810,8255,1678,1015,5304,7603,2023,3855,1969,2103,4690,2072,8814,1322,3915,7015,1353,4379,560,5483,8254,9322,216,6506,1090,4444,5098,7510,8417,6169,4374,426,6153,595,9061,3346,1594,6515,2439,2792,3582,5612,9964,7100,7309,4428,3787,3526,7269,10347,6764,4657,8065,9751,1634,547,4263,6976,9817,5128,6772,1855,3000,9,8748,4273,8011,2454,4130,9111,2331,7972,8058,8414,1659,2505,5773,6113,3166,178,4369,3609,2568,2681,601,1986,7071,10073,9520,7506,7308,2303,4606,8334,494,5290,4245,9195,8379,8052,1625,4147,7937,5438,1683,4216,8002,2694,1574,9143,8758,5239,2657,2252,9437,6689,5638,2625,7512,770,3813,10360,3599,123,3931,905,2549,9370,5757,10439,3693,6653,10259,2680,6444,5564,8587,5649,1861,9998,3179,2812,7968,9387,5524,4782,8481,3566,7355,4349,5296,9281,2500,1374,8094,9893,9667,6121,4347,3203,3414,842,10128,5965,2061,5601,2623,6374,5133,7766,8950,2205,7813,4030,6287,7995,5046,1240,7131,8163,593,1391,5196,8351,9940,5430,2877,9668,1553,5469,9587,3355,3353,10176,2410,2396,9160,3549,7307,7418,9091,6078,5376,802,7393,3617,10465,2693,2732,2482,2848,9333,4125,7616,6765,7800,7033,10156,5195,8288,4752,6488,667,793,1741,1575,6105,5722,8902,10342,3387,7171,4284,474,9512,9734,3548,3097,281,6964,3141,210,8796,7606,98,3370,8773,1348,2108,849,5288,8708,6762,9217,7036,9293,5504,7194,467,3275,5531,5878,4892,4089,2768,8503,2075,6532,3263,1935,5721,160,9580,2003,5783,8015,2412,6148,2705,2981,5220,889,6248,2485,4082,7817,8014,1050,2523,10110,3547,10310,8668,9664,9581,2175,3208,6795,9954,5788,242,5054,382,3364,116,4183,5271,2138,1825,2305,8487,5300,3532,3979,806,1957,5112,5755,2824,4560,2791,7588,4501,10107,8046,4255,988,7946,7756,351,833,7599,5073,1317,3751,7522,7666,1783,8489,375,4603,791,5699,5656,9631,9750,10,1812,8507,9761,2263,3505,7605,3292,2594,1795,755,8151,1410,3196,6538,7058,2167,5899,423,3753,6260,2487,2522,5328,5506,9907,6047,10237,6168,2664,2215,1657,2917,5121,5623,6736,436,2124,9625,6543,6004,6877,3891,10076,10384,8704,6639,8611,9104,5301,2306,1511,6453,8871,9950,6780,1185,6769,164,3990,7710,4354,9268,8170,1405,490,4045,4944,4926,8608,7098,10155,6282,7939,4391,4211,1036,3954,4408,7920,963,1154,8724,7049,838,4539,8053,8018,10477,2136,2554,9681,6781,10263,6272,6451,2921,6708,9695,9203,989,4756,6729,794,7531,596,403,440,8017,5425,2054,2561,4848,5323,43,4914 tcga_RSEM_gene_fpkm > BRCA_tcga_RSEM_gene_fpkm

cp BRCA_tcga_RSEM_gene_fpkm /home/erpl/Project/5.Classification_of_TCGA_samples/Custom
cp BRCA_tcga_RSEM_gene_fpkm /home/erpl/Project/5.Classification_of_TCGA_samples/msigdb

cd /home/erpl/Project/5.Classification_of_TCGA_samples/Custom
#merging custom set with ensembl id
#ensembl_ID.txt contains IDs for all human genes and their transcripts downloaded from Ensembl Biomart
R
a=read.table(file="custom_sig_after_cleanup_34gene.txt",header=F,sep="\t")
b=read.table(file="ensembl_ID.txt",header=T,sep="\t")
c=merge(a,b,by.x="V1",by.y="Gene_name")
write.table(c,file="custom_m_ID",append=F,sep="\t",row.names=F,quote=F)
q()
n
cut -f 1,2 custom_m_ID|sort|uniq|grep -v "V1">custom_signature_34_genes.txt

#subset for custom genes expression from whole set

for i in `cat custom_signature_34_genes.txt|cut -f2`
do
grep "$i" BRCA_tcga_RSEM_gene_fpkm >> BRCA_fpkm_custom
done

head -1 BRCA_tcga_RSEM_gene_fpkm>header
cat header BRCA_fpkm_custom>BRCA_tcga_RSEM_gene_custom_signature_34genes

datamash transpose < BRCA_tcga_RSEM_gene_custom_signature_34genes>BRCA_tcga_RSEM_gene_custom_signature_34genes_trans
R
a=read.table(file="BRCA_tcga_RSEM_gene_custom_signature_34genes_trans",header=T,sep="\t")
sum_a=cbind(a,total=rowSums(a[,-1]))
write.table(sum_a,file="sum_custom",append=F,sep="\t",row.names=F,quote=F)
q()
n
sort -k 36 -g -r sum_custom>sum_custom_sorted

#For obtaining top and bottom 25% samples
sed '277,825d' sum_custom_sorted > sum_custom_top_bottom_25
cut -f 1-35 sum_custom_top_bottom_25>sum_top_bottom_custom25_cluster
grep -v "sample" sum_top_bottom_custom25_cluster|cut -f 1 |head -275>hypoxia_cluster_custom
grep -v "sample" sum_top_bottom_custom25_cluster|cut -f 1 |tail -275>normoxia_cluster_custom

#For obtaining top and bottom 20% samples
sed '222,880d' sum_custom_sorted > sum_custom_top_bottom_20
cut -f 1-35 sum_custom_top_bottom_20>sum_top_bottom_custom20_cluster
grep -v "sample" sum_top_bottom_custom20_cluster|cut -f 1 |head -220>hypoxia_cluster_custom20
grep -v "sample" sum_top_bottom_custom20_cluster|cut -f 1 |tail -220>normoxia_cluster_custom20

#For obtaining top and bottom 10% samples
sed '112,990d' sum_custom_sorted > sum_custom_top_bottom_10
cut -f 1-35 sum_custom_top_bottom_10>sum_top_bottom_custom10_cluster
grep -v "sample" sum_top_bottom_custom10_cluster|cut -f 1 |head -110>hypoxia_cluster_custom10
grep -v "sample" sum_top_bottom_custom10_cluster|cut -f 1 |tail -110>normoxia_cluster_custom10


#classification of msigdb
cd /home/erpl/Project/5.Classification_of_TCGA_samples/msigdb

#merging msigdb set with ensembl id
R
a=read.table(file="msigdb_signature.txt",header=F,sep="\t")
b=read.table(file="ensembl_ID.txt",header=T,sep="\t")
c=merge(a,b,by.x="V1",by.y="Gene_name")
write.table(c,file="hypoxia_hallmark_msigdb.txt",append=F,sep="\t",row.names=F,quote=F)
q()
n

cut -f 1,2 hypoxia_hallmark_msigdb.txt|sort|uniq|grep -v "V1">msigdb_signature_ensembl_ID.txt


for i in `cat msigdb_signature_ensembl_ID.txt|cut -f1`
do
grep "$i" BRCA_tcga_RSEM_gene_fpkm >> BRCA_fpkm_msigdb
done

head -1 BRCA_tcga_RSEM_gene_fpkm>header
cat header BRCA_fpkm_msigdb>BRCA_tcga_RSEM_gene_hypoxia_hallmark_msigdb


datamash transpose < BRCA_tcga_RSEM_gene_hypoxia_hallmark_msigdb>BRCA_tcga_RSEM_gene_hypoxia_hallmark_msigdb_trans

R
b=read.table(file="BRCA_tcga_RSEM_gene_hypoxia_hallmark_msigdb_trans",header=T,sep="\t")
sum_b=cbind(b,total=rowSums(b[,-1]))
write.table(sum_b,file="sum_msigdb",append=F,sep="\t",row.names=F,quote=F)
q()
sort -k 202 -r sum_msigdb>sum_msigdb_sorted

#For obtaining top and bottom 25% samples
sed '277,825d' sum_msigdb_sorted > sum_msigdb_top_bottom_25
cut -f 1-201 sum_msigdb_top_bottom_25>sum_top_bottom_msigdb25_cluster
grep -v "sample" sum_msigdb_top_bottom_25_cluster|cut -f 1 |head -275>hypoxia_cluster_msigdb
grep -v "sample" sum_msigdb_top_bottom_25_cluster|cut -f 1 |tail -275>normoxia_cluster_msigdb

#For obtaining top and bottom 20% samples
sed '222,880d' sum_msigdb_sorted > sum_msigdb_top_bottom_20
cut -f 1-201 sum_msigdb_top_bottom_20>sum_top_bottom_msigdb20_cluster
grep -v "sample" sum_msigdb_top_bottom_20_cluster|cut -f 1 |head -220>hypoxia_cluster_msigdb20
grep -v "sample" sum_msigdb_top_bottom_20_cluster|cut -f 1 |tail -220>normoxia_cluster_msigdb20

#For obtaining top and bottom 10% samples
sed '112,990d' sum_msigdb_sorted > sum_msigdb_top_bottom_10
cut -f 1-201 sum_msigdb_top_bottom_10>sum_top_bottom_msigdb10_cluster
grep -v "sample" sum_msigdb_top_bottom_10_cluster|cut -f 1 |head -110>hypoxia_cluster_msigdb10
grep -v "sample" sum_msigdb_top_bottom_10_cluster|cut -f 1 |tail -110>normoxia_cluster_msigdb10
###################################################################################################################################################################

#Code for kmeans clustering and ROC

for i in `ls -1|grep -v "kmeans.r"| cut -f 4 -d '_'`; do Rscript kmeans.r $i; done

#Rscript kmeans.r
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#S=read.table(file=paste("sum_top_bottom_",args[1],"_cluster",sep=""),header=TRUE,sep="\t")
#prcomp(S[,-1], center = TRUE,scale. = TRUE)->PCAS
#data.frame(PC1=PCAS$x[,1],PC2=PCAS$x[,2],sam=S$sample)->hun
#rownames(S)<-S$sample
#kmeans(S[,2:ncol(S)], 2, nstart = 20)->tSclusters
#as.data.frame(tSclusters$cluster)->Q
#rownames(Q)->Q$sample
#colnames(Q)<-c("Cluster","sample")
#merge(hun,Q,by.x="sam",by.y="sample")->hunQ
#tiff(paste("kmeans_",args[1],".tiff",sep=""),units="in",height=5,width=5,res=300)
#plot(hunQ$PC1,hunQ$PC2,col=ifelse(hunQ$Cluster==2,"red","blue"),pch=16,xlab="PC1",ylab="PC2")
#dev.off()
#nrow(S)/2->n
#data.frame(Sample=S$sample,Truth=c(rep(0,n),rep(1,n)))->ROC
#merge(ROC,hunQ,by.x="Sample",by.y="sam")->ROCp
#write.table(file=paste("truth_predicted_ROC_",args[1],".txt",sep=""),ROCp,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

read.table(file="truth_predicted_ROC.txt",header=FALSE)->T
T$V5-1->T$V6
library(ROCR)
pred <- prediction(T$V6, T$V2)
perf <- performance(pred,"tpr","fpr")
tiff("",units="in",height=5,width=5,res=300)
plot(perf,colorize=TRUE)
dev.off()
q()
n

#For calculating classification evaluation parameters
cut -f 2 truth_predicted_ROC|sed 's/0/abnormal/g'|sed 's/1/normal/g' > Truth.txt
cut -f 5 truth_predicted_ROC|sed 's/2/abnormal/g'|sed 's/1/normal/g' > Predict.txt
R
read.table(file="Truth.txt",header=FALSE,as.is=TRUE)->T
read.table(file="Predict.txt",header=FALSE,as.is=TRUE)->P
xtab <- table(P$V1, T$V1)
# load Caret package for computing Confusion matrix
library(caret)
confusionMatrix(xtab)
q()
n




######################################################################################################################################################################
############################################6. Differential Exon Expression Analysis####################################################################################
######################################################################################################################################################################


#For calculating Isoform Fractions

#Splitting the custom clusters into 5 parts
split -l 25 normoxia_cluster_custom10 NormC_custom.
split -l 25 hypoxia_cluster_custom10 HypC_custom.

#Creating abundance.tsv for every breast cancer tumor sample
for i in `cat HypC_custom.aa`
do
echo $i
est=`cat datafile | tr "\t" "\n"|awk -v tar=$i '$1==tar{print NR}'`
cut -f $est tcga_Kallisto_est_counts > temp."$i".counts
cut -f $est tcga_Kallisto_tpm > temp."$i".tpm
paste trans.list temp."$i".counts temp."$i".tpm > temp."$i".kall
mkdir -p estimates/"$2"_"$i"
Rscript WriteKalAb.r temp."$i".kall "$2"_"$i"
rm temp."$i".counts temp."$i".tpm temp."$i".kall
done

#WriteKalAb.r Rscipt
#args = commandArgs(trailingOnly=TRUE)
#read.table(file="kall_format",header=FALSE)->E
#read.table(file=args[1],header=TRUE,as.is=TRUE)->M
#merge(E,M,by.y="sample",by.x="V1")->MEI
#data.frame(target_id=MEI$V1,length=MEI$V2,eff_length=MEI$V3,est_counts=(2^MEI[,4])-1,tpm=(2^MEI[,5])-0.001)->X
#X$tpm[X$tpm<0]<-rep(0,length(X$tpm[X$tpm<0]))
#write.table(file=paste("estimates/",args[2],"/abundance.tsv",sep=""),X,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


#Same was repeated for all other split clusters

#Calculating Isoform Fraction for every transcript
for j in `ls -1 estimates/|cut -f 2 -d '/'`
do
sed 's/\./\t/' estimates/"$j"/abundance.tsv|grep -v "target_id" > "$j"_abund
Rscript getIF.r "$j"_abund
done

#getIF.r Rscript 
args = commandArgs(trailingOnly=TRUE)
read.table(file="ensembl_ID.txt",header=TRUE,as.is=TRUE)->M
read.table(file=args[1],header=FALSE,as.is=TRUE,fill=TRUE)->E
merge(M,E,by.x="Transcript_stable_ID",by.y="V1")->ME
as.data.frame(aggregate(ME$V5,list(ME$Gene_stable_ID),sum))->N
merge(ME,N,by.x="Gene_stable_ID",by.y="Group.1")->MEN
MEN$V5/MEN$x->MEN$IF
data.frame(T=MEN$Transcript_stable_ID,IF=MEN$IF)->X
write.table(file=paste(args[1],".IF",sep=""),X,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

paste HypC_*_abund.IF|cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220 > IF_HypC
paste NormC_*_abund.IF|cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220 > IF_NormC

#Getting mean of isoform fraction and difference in isoform fraction dIF
read.table(file="IF_HypC",header=FALSE)->H
apply(H[,2:111], 1, mean)->HX
read.table(file="IF_NormC",header=FALSE)->N
apply(N[,2:111], 1, mean)->NX
HX-NX->DIF
na.omit(data.frame(Trans=H$V1,DIF=DIF))->Q
#as.numeric(abs(quantile(Q$DIF,0.01,na.rm=TRUE)))
#take the FDR<0.01 cut-off
unique(na.omit(Q$Trans[abs(Q$DIF)>as.numeric(abs(quantile(Q$DIF,0.01,na.rm=TRUE)))]))->UTrans
read.table(file="ensembl_ID.txt",header=TRUE,as.is=TRUE)->E
data.frame(table(E$Gene_stable_ID[E$Transcript_stable_ID %in% UTrans]))->G2
G2[G2$Freq>1,]->G3
write.table(file="dIF_C_genes_list.txt",G3$Var1,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

#creating file with gene_ID, transcript_ID and dIF

Q[as.character(Q$Trans) %in% as.character(UTrans),]->Q2
merge(Q2,E,by.x="Trans",by.y="Transcript_stable_ID")->Q2E
Q2E[Q2E$Gene_stable_ID %in% G3$Var1,]->Q3
write.table(file="dIF_C_genes_list_IF.txt",Q3,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

#all the above steps were repeated for msigdb set

#for plotting boxplot for candidate genes
I1<-"ENST00000260947"
I2<-"ENST00000619009"

H[H$V1==I1,]->HI1
H[H$V1==I2,]->HI2

N[N$V1==I1,]->NI1
N[N$V1==I2,]->NI2

data.frame(Factor=c(rep("B_HI1",length(HI1)-1),rep("D_HI2",length(HI2)-1),rep("A_NI1",length(NI1)-1),rep("C_NI2",length(NI2)-1)),Expr=c(as.numeric(HI1[1,2:111]),as.numeric(HI2[1,2:111]),as.numeric(NI1[1,2:111]),as.numeric(NI2[1,2:111])))->B

boxplot(B$Expr~B$Factor,noth=TRUE,col=c("red","blue","red","blue"),ylab="Percent Isoform Usage",names=c("Normoxia","Hypoxia","Normoxia","Hypoxia"))



####################################################################################################################################################################
######################################################7. Correlation Test For Methylation#############################################################################
####################################################################################################################################################################

cat HiSeqV2_exon|grep -v "Sample"|cut -f 1|awk '{print $0,$0}'|sed 's/:/\t/'|sed 's/-/\t/'|sed 's/:/\t/'|grep "chr"|sed 's/chr//'|sort -k1n,1 -k2n,2 > pre_subtract_exons.bed

#for extracting transcripts that have at least 2 unique transcripts
R
unqlen <- function(temp_F) {
  length(unique(temp_F))
}

read.table(file="biomart_exons.bed",header=FALSE)->E
aggregate(E$V5,list(E$V6),unqlen)->EU
EU$Group.1[EU$x>1]->EUT
E[E$V6 %in% EUT,]->ET
write.table(file="new_exons.bed",ET,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
q()
n

#for obtaining those genes which have at least one exon which is getting differentially spliced between the 2 transcripts 
for gene in `cut -f 6 new_exons.bed|sort -u`
do
count=1
for trans in `grep "$gene" new_exons.bed|cut -f 5|sort -u`
do
echo $gene $trans
grep "$trans" new_exons.bed > trans."$count".bed
count=`echo $count|awk '{print $1+1}'`
done
tran1=`cat trans.1.bed|wc -l`
tran2=`cat trans.2.bed|wc -l`
if [ $tran1 -gt 1 ]
then
if [ $tran2 -gt 1 ]
then
bedtools subtract -A -a trans.1.bed -b trans.2.bed >> subtracted_regions.bed
bedtools subtract -A -a trans.2.bed -b trans.1.bed >> subtracted_regions.bed
fi
fi
done


bedtools intersect -f 1 -F 1 -a pre_subtract_exons.bed -b subtracted_regions.bed -wa|sort -u > exons.bed

for i in `ls -1 dIF_C.*`
do
nohup ./getCors.sh $i &
done

#for msigdb

grep -v "Gene" diff_spliced_genes_msigdb.txt|awk '{print $4,$2,$3,$1}'|sed 's/ /\t/g'>diff_spliced_genes_msigdb.bed
split -l 50 diff_spliced_genes_msigdb.bed dIF_msigdb.*
for i in `ls -1 dIF_msigdb.*`
do
nohup ./getCors.sh $i &
done  


#getCors.sh script
#for gene in `cat $1|sed 's/\t/_/g'`
#do
#echo $gene
#echo $gene|sed 's/_/\t/g' > "$gene".bed
#bedtools intersect -a "$gene".bed -b illuminaMethyl450_hg19_GPL16304_TCGAlegacy.bed -wb|cut -f 5-8|sort -k1n,1 -k2n,2 > probes_in_"$gene".bed
#bedtools intersect -a "$gene".bed -b exons.bed -wb|cut -f 5-8|sort -k1n,1 -k2n,2 > exons_in_"$gene".bed
#for pair in `bedtools closest -k 5 -d -a exons_in_"$gene".bed -b probes_in_"$gene".bed -wa -wb|awk '$10<1000{print $0}'|cut -f 4,8,9|sed 's/\t/_/g'|cut -f 2 -d ' '`
#do
#exonP=`echo $pair|cut -f 1 -d '_'`
#methP=`echo $pair|cut -f 2 -d '_'`
#head -1 HiSeqV2_exon > exon."$gene".temp
#grep "$exonP" HiSeqV2_exon >> exon."$gene".temp
#head -1 meth_BRCA_PANCAN.txt > temp."$gene".meth
#grep "$methP" meth_BRCA_PANCAN.txt >> temp."$gene".meth
#Rscript getCor.r $gene
#rm temp."$gene".meth exon."$gene".temp
#done
#rm "$gene".bed exons_in_"$gene".bed probes_in_"$gene".bed 
#mv "$gene".Cor_out.txt cors_folder
#done

#getCor.r script
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#read.table(file=paste("exon.",args[1],".temp",sep=""),header=FALSE)->E
#read.table(file=paste("temp.",args[1],".meth",sep=""),header=FALSE,fill=TRUE)->M
#data.frame(t(E)[-1,])->tE
#data.frame(t(M)[-1,])->tM
#merge(tE,tM,by.x="X1",by.y="X1")->tEM
##plot(as.numeric(as.character(tEM$X2.x)),as.numeric(as.character(tEM$X2.y)),xlab="Methylation",ylab="Expression",pch=16)
#cor.test(as.numeric(as.character(tEM$X2.x)),as.numeric(as.character(tEM$X2.y)))->P
#cor.test(as.numeric(as.character(tEM$X2.x)),as.numeric(as.character(tEM$X2.y)),method="kendall")->K
#data.frame(P$estimate,P$p.value,K$estimate,K$p.value,t(E)[1,2],t(M)[1,2],args[1])->X
#write.table(file=paste(args[1],".Cor_out.txt",sep=""),X,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

####################################################################################################################################################################
#########################################################Survival Analysis##############################################################################################
####################################################################################################################################################################

mkdir Survival_analysis
cd Survival_analysis
cp /home/erpl/Project/4.Collection_of_TCGA_data/Survival_SupplementalTable_S1_20171025_xena_sp ../Survival_analysis
cp /home/erpl/Project/5.Classification_of_TCGA_samples/tumor_sample_IDs ../Survival_analysis/
for i in `cat tumor_sample_IDs`; do  grep "$i" Survival_SupplementalTable_S1_20171025_xena_sp ; done>BRCA_survival_data
head -1 Survival_SupplementalTable_S1_20171025_xena_sp > header
cat header BRCA_survival_data > BRCA_survival_data_withHeader
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/hypoxia_cluster_custom10 ../Survival_analysis/
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/normoxia_cluster_custom10 ../Survival_analysis/
awk '{$(NF+1)=1;}1' hypoxia_cluster_custom10 > hypoxia_cluster_custom10_1
awk '{$(NF+1)=2;}1' normoxia_cluster_custom10 >normoxia_cluster_custom10_2
cat hypoxia_cluster_custom10_1 normoxia_cluster_custom10_2 > combined_clusters_custom
sed -i 's/ /\t/g' combined_clusters_custom
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_custom",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom_BRCA_survival",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 custom_BRCA_survival |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table
R
library(survival)
library(survminer)
a=read.table(file="surv_table",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)
tiff("survival_curve_custom.tiff",units="in",height=10,width=10,res=300)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Hypoxia", "Normoxia"), legend.title="Condition",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.15)
dev.off()
fit.coxph <- coxph(Surv(time_alive_or_dead, vital_status)~Condition, 
                   data = surv)
tiff("HR_plot_custom.tiff",units="in",height=6,width=6,res=300)
ggforest(fit.coxph, data = surv)
dev.off()
q()
n

#Survival curve for msigdb clusters
cp /home/erpl/Project/5.Classification_of_TCGA_samples/msigdb/normoxia_cluster_msigdb10 ../Survival_analysis/
cp /home/erpl/Project/5.Classification_of_TCGA_samples/msigdb/hypoxia_cluster_msigdb10 ../Survival_analysis/
awk '{$(NF+1)=1;}1' hypoxia_cluster_msigdb10 > hypoxia_cluster_msigdb10_1
awk '{$(NF+1)=2;}1' normoxia_cluster_msigdb10 >normoxia_cluster_msigdb10_2
cat hypoxia_cluster_msigdb10_1 normoxia_cluster_msigdb10_2 > combined_clusters_msigdb
sed -i 's/ /\t/g' combined_clusters_msigdb
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_msigdb",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="msigdb_BRCA_survival",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 msigdb_BRCA_survival |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_msigdb
R
a=read.table(file="surv_table_msigdb",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)
tiff("surv_plot_msigdb.tiff",units="in",height=10,width=10,res=300)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Hypoxia", "Normoxia"), legend.title="Condition",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.15)
dev.off()
fit.coxph <- coxph(Surv(time_alive_or_dead, vital_status)~Condition, 
                   data = surv)
tiff("HR_plot_msigdb.tiff",units="in",height=6,width=6,res=300)
ggforest(fit.coxph, data = surv)
dev.off()
q()
n

#Top 5% of samples for survival plot
#For obtaining top and bottom 5% samples
#For obtaining top and bottom 10% samples
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/sum_custom_sorted
grep -v "sample" sum_custom_sorted|sed '56,1044d'> sum_custom_top_bottom_5
cut -f 1-35 sum_custom_top_bottom_5>sum_top_bottom_custom5_cluster
cut -f 1 sum_top_bottom_custom5_cluster|head -55>hypoxia_cluster_custom5
cut -f 1 sum_top_bottom_custom5_cluster|tail -55>normoxia_cluster_custom5

awk '{$(NF+1)=1;}1' hypoxia_cluster_custom5 > hypoxia_cluster_custom5_1
awk '{$(NF+1)=2;}1' normoxia_cluster_custom5 >normoxia_cluster_custom5_2
cat hypoxia_cluster_custom5_1 normoxia_cluster_custom5_2 > combined_clusters_custom5
sed -i 's/ /\t/g' combined_clusters_custom5

R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_custom5",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom5_BRCA_survival",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 custom5_BRCA_survival |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_5
R
library(survival)
a=read.table(file="surv_table_5",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)

#Top 1% of samples for survival plot
grep -v "sample" sum_custom_sorted|sed '12,1088d'> sum_custom_top_bottom_1
cut -f 1-35 sum_custom_top_bottom_1>sum_top_bottom_custom1_cluster
cut -f 1 sum_top_bottom_custom1_cluster|head -11>hypoxia_cluster_custom1
cut -f 1 sum_top_bottom_custom1_cluster|tail -11>normoxia_cluster_custom1
awk '{$(NF+1)=1;}1' hypoxia_cluster_custom1 > hypoxia_cluster_custom1_1
awk '{$(NF+1)=2;}1' normoxia_cluster_custom1 >normoxia_cluster_custom1_2
cat hypoxia_cluster_custom1_1 normoxia_cluster_custom1_2 > combined_clusters_custom1
sed -i 's/ /\t/g' combined_clusters_custom1

R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_custom1",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom1_BRCA_survival",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 custom1_BRCA_survival |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_1
R
library(survival)
a=read.table(file="surv_table_1",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)

##Survival curve according to gene expression
datamash transpose < sum_custom_sorted_withHeader > sum_custom_sorted_withHeader_trans
for i in `cat sum_custom_sorted_withHeader_trans`
do
cut -f 2 $i>gene
echo gene
done

##For alkbh5
cut -f 1,2 sum_custom_sorted_withHeader|sort -k 2 -r|grep -v "sample"|cut -f 1|sed '111,989d'>top_bottom_10_ALKBH5
head -110 top_bottom_10_ALKBH5>hypoxia_ALKBH5_10
tail -110 top_bottom_10_ALKBH5>normoxia_ALKBH5_10


awk '{$(NF+1)=1;}1' hypoxia_ALKBH5_10 > hypoxia_ALKBH5_10_1
awk '{$(NF+1)=2;}1' normoxia_ALKBH5_10 >normoxia_ALKBH5_10_2
cat hypoxia_ALKBH5_10_1 normoxia_ALKBH5_10_2 > combined_clusters_alkbh5
sed -i 's/ /\t/g' combined_clusters_alkbh5
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_alkbh5",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="alkbh5_BRCA_survival",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 alkbh5_BRCA_survival |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_alkbh5
R
a=read.table(file="surv_table_alkbh5",header=T,sep="\t")
surv=as.data.frame(a)
library(survival)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)
tiff("surv_plot_msigdb.tiff",units="in",height=10,width=10,res=300)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Hypoxia", "Normoxia"), legend.title="Condition",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for Breast Cancer Survival", 
           risk.table.height=.15)
dev.off()


#25% alkbh5
cut -f 1,2 sum_custom_sorted_withHeader|sort -k 2 -r|grep -v "sample"|cut -f 1|sed '276,824d'>top_bottom_25_ALKBH5
head -275 top_bottom_25_ALKBH5>hypoxia_ALKBH5_25
tail -275 top_bottom_25_ALKBH5>normoxia_ALKBH5_25

awk '{$(NF+1)=1;}1' hypoxia_ALKBH5_25 > hypoxia_ALKBH5_25_1
awk '{$(NF+1)=2;}1' normoxia_ALKBH5_25 >normoxia_ALKBH5_25_2
cat hypoxia_ALKBH5_25_1 normoxia_ALKBH5_25_2 > combined_clusters_alkbh5_25
sed -i 's/ /\t/g' combined_clusters_alkbh5_25
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_alkbh5_25",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="alkbh5_BRCA_survival_25",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 alkbh5_BRCA_survival_25 |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_alkbh5_25
R
a=read.table(file="surv_table_alkbh5_25",header=T,sep="\t")
surv=as.data.frame(a)
library(survival)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)

#25% custom

cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/hypoxia_cluster_custom ../Survival_analysis/
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/normoxia_cluster_custom ../Survival_analysis/
awk '{$(NF+1)=1;}1' hypoxia_cluster_custom > hypoxia_cluster_custom25_1
awk '{$(NF+1)=2;}1' normoxia_cluster_custom >normoxia_cluster_custom25_2
cat hypoxia_cluster_custom25_1 normoxia_cluster_custom25_2 > combined_clusters_custom25
sed -i 's/ /\t/g' combined_clusters_custom25
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_custom25",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom_BRCA_survival_25",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 custom_BRCA_survival_25 |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_custom25
R
library(survival)
a=read.table(file="surv_table_custom25",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)

#20% custom
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/hypoxia_cluster_custom20 ../Survival_analysis/
cp /home/erpl/Project/5.Classification_of_TCGA_samples/Custom/normoxia_cluster_custom20 ../Survival_analysis/
awk '{$(NF+1)=1;}1' hypoxia_cluster_custom20 > hypoxia_cluster_custom20_1
awk '{$(NF+1)=2;}1' normoxia_cluster_custom20 >normoxia_cluster_custom20_2
cat hypoxia_cluster_custom20_1 normoxia_cluster_custom20_2 > combined_clusters_custom20
sed -i 's/ /\t/g' combined_clusters_custom20
R
a=read.table(file="BRCA_survival_data_withHeader",header=T,sep="\t")
b=read.table(file="combined_clusters_custom20",header=F,sep="\t")
c=merge(a,b,by.x="sample",by.y="V1")
write.table(c,file="custom_BRCA_survival_20",sep="\t",append=F,row.names=F,quote=F)
q()
n
cut -f 4,5,14,16,17,35 custom_BRCA_survival_20 |sed 's/NA/0/g'|sed 's/Alive/1/g'|sed 's/Dead/2/'|sed 's/FEMALE/2/g'|sed 's/MALE/1/g'|awk '{print $1,$2,$3,$4+$5,$6}'|sed 's/ /\t/g'| sed '1s/0/time_alive_or_dead/g'|sed 's/V2/Condition/'>surv_table_custom20
R
library(survival)
a=read.table(file="surv_table_custom20",header=T,sep="\t")
surv=as.data.frame(a)
sfit <- survfit(Surv(time_alive_or_dead, vital_status)~Condition, data=surv)










