"""
Global variables called during the GORi data preparation step.

@ May 2021
@ Asloudj Yanis
"""

from _scripts.misc import read_csv

# set of genes available for testing:
test_genes1 = "LY6E,IFIT1,OAS1,IFIT1,IFIT3,OAS3,IFIT3,OAS1,OASL,LOC129607,ISG15,HERC5,OAS1,MX1,BATF2,LAMP3,IFI44L,XAF1,OASL,IFI44,OAS2,TRIM6,HES4,OTOF,FLJ20035,IFITM3,IFIT3,CXCL10,EPSTI1,SERPING1,LOC26010,OAS2,RSAD2,RTP4,"
test_genes1 = test_genes1.split(",")
test_genes2 = "NCKAP1L,THEMIS,HLA-DPB1,JUN,NFKBIA,PTPN6,TLR5,HLA-DRA,CD4,PRKCQ,PLCG1,CD247,TRBC1,FCGR1B,CCR7,TLR6,TLR10,FYB,TLR2,FCER1G,TRAT1,LILRB2,LILRB1,CLEC7A,BTLA,FYN,FGR,CD86,GRAP2,HLA-DQB1,ZAP70,PAG1,MNDA,IRAK3,LCK,DUSP6,PTPN22,PAK1,CTSH,CTLA4,LCP2,TLR4,CSK,FOS,LAT2,CARD11,HLA-DRB1,CHUK,CD24,UBASH3A,MAP2K6,TLR7,LY96,HCK,VASP,LYN,LIME1,MAPK14,MYD88,INPP5D,ITK,TRAC,SLA2,KLRK1,HLA-DMB,RPS6KA1,HLA-DQA1,MEF2C,SYK,MEF2A,TLR1,PTPRC,HLA-DOA,HLA-DOB,HLA-DMA,HLA-DPA1,CD19,CD3D,CD3E,BTK,CD3G"
test_genes2 = test_genes2.split(",")
# test_genes3 = read_csv("./data/GSE31684_Intersect_genes.csv")
pig_test_genes = ['A0A287A243', 'A0A286ZP59', 'A0A287AXZ5', 'I3LS66', 'F1RKA6', 'F1S8Q0', 'A0A286ZYB0', 'I3LQQ9', 'F1RVZ4', 'A0A287A147', 'F1RJC9', 'K7GNR7', 'A0A287A1X9', 'A0A286ZLE1', 'I3LNM4', 'I3LHQ8', 'A0A5G2R0C6', 'A0A287B828', 'F1S2D4', 'A0A480MHE7', 'P51779', 'F1RLC1', 'A0A286ZWY1', 'A0A287BC48', 'I3LBM6', 'A0A5G2R8A1', 'A0A287B934', 'A0A287BPX6', 'F1S866', 'P60899', 'I3LMK3', 'A0A287AL63', 'A0A287AA05', 'A0A5G2R2L7', 'I3LQH7', 'A0A5G2QTP9', 'A0A5G2QIZ6', 'F1S5C4', 'A0A5G2QU22', 'A0A5G2QVC6', 'I3L6G1', 'A0A287AC75', 'F1S6H0', 'A0A5G2QUJ5', 'A0A287BCU5', 'I3LBF9', 'F1SM17', 'I3LCQ9', 'A0A286ZWZ8', 'F1RHR8', 'Q6TYI6', 'F1RZ33', 'A0A287B5G8', 'F1RQB7', 'A0A5G2QJQ3', 'F1SMJ1', 'F1SIU1', 'F1SRH7', 'A0A287AJT7', 'A0A287AUN9', 'I3LJD2', 'A0A286ZLS4', 'I3LA25', 'E1U8C4', 'A0A287BHP2', 'A0A5G2RKA4', 'F1S8S0', 'A0A287A3D1', 'F1RQU1', 'K7GPB4', 'A0A287AFR9', 'F1RU77', 'A0A287A316', 'A0A5G2QH83', 'F1RVX6', 'F1RF12', 'F1SAM0', 'I3LJG4', 'A0A286ZK65', 'I3L8U7', 'A0A287AI58', 'A0A5G2QZV9', 'I3LRY2', 'A0A287BHC5', 'A0A286ZM79', 'A0A5G2QSH2', 'F1SAH5', 'I3LRY4', 'A0A5G2R802', 'Q52NJ1', 'A0A287BMZ8', 'F1RXV9', 'A0A5G2QY32', 'Q6YT47', 'F1SCF3', 'A0A5G2R2J2', 'A0A286ZU26', 'A0A5G2R181', 'F1S544', 'F1S8J8', 'F1SKG1', 'A0A286ZQV5', 'A0A5G2QIS7', 'I3LG77', 'A0A287ALH4', 'F1S6W6', 'A0A286ZPJ7', 'A0A5G2QSK0', 'F1SDX0', 'A0A287A8L1', 'A0A287AQT2', 'K7GN73', 'A0A287A8Q7', 'F1RG22', 'F1RQG9', 'F1RRK8', 'A0A288CFY2', 'A0A286ZHZ3', 'A0A5G2QY29', 'A0A287AD64', 'F1SKU0', 'A0A5G2QCN4', 'A0A287AR65', 'Q02038', 'F1RXR3', 'F1S8C9', 'F1RSK1', 'A0A5G2QEE3', 'F1S397', 'I3LGV0', 'K7GKA9', 'A0A287A042', 'A0A481CED7', 'I3LJS5', 'B0M1M6', 'A0A287BB81', 'A0A286ZL51', 'A0A287A1L0', 'A0A287A8Y8', 'F1SEK6']

# the GO terms aspect of interest:
target_aspect = ["biological_process"]
# the genes of interest:
target_genes = test_genes2
# the ontologies of interest:
target_onto = ["GO", "R-", "HP"]
# genes are symbols instead of UniProtKB IDS:
symbol = True
# If the GAF file is downloaded from the Gene Ontology database:
download = True # more details below #

## relative pathes on your computer:
data_path = "./data"
results_path = "./results"

## URLs to download data files from and their respective resulting files:
# GAF file is downloaded automatically by a function call.
"""
The GAF file of interest:
    1. Select a GAF file from http://current.geneontology.org/annotations/index.html 
    2. Copy-paste it:
"""

# Indicate the relative path leading to the GAF file of interest:
# Download a GAF file from GO:
gaf_name = "goa_human.gaf.gz" # http://current.geneontology.org/annotations/index.html
gaf_file = data_path + "/" + gaf_name[:-3]
# Or use a local GAF file :
# gaf_file = "my_documents/my_gaf_file.gaf"

go_obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
go_obo_file = "%s/go-basic.obo" % data_path

reactome_hierarchy_url = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
reactome_hierarchy_file = "%s/reactome_hierarchy.txt" % data_path

reactome_label_url = "https://reactome.org/download/current/ReactomePathways.txt"
reactome_label_file = "%s/reactome_label.txt" % data_path

reactome_annotation_url = "https://reactome.org/download/current/UniProt2Reactome.txt"
reactome_annotation_file = "%s/reactome_annotation.txt" % data_path

# reactome obo file is not downloaded but generated.
reactome_obo_file = "%s/reactome.obo" % data_path

hpo_annotation_url = "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
hpo_annotation_file = "%s/hpo_annotation.txt" % data_path

hpo_obo_url = "http://purl.obolibrary.org/obo/hp.obo"
hpo_obo_file = "%s/hp.obo" % data_path

# results files
items_file = "%s/items.csv" % results_path
rules_file = "%s/rules.csv" % results_path
