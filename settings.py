"""
Global variables called during the GORAFI data preparation step.

@ May 2021
@ Asloudj Yanis
"""

# set of genes available for testing:
test_genes1 = "LY6E,IFIT1,OAS1,IFIT1,IFIT3,OAS3,IFIT3,OAS1,OASL,LOC129607,ISG15,HERC5,OAS1,MX1,BATF2,LAMP3,IFI44L,XAF1,OASL,IFI44,OAS2,TRIM6,HES4,OTOF,FLJ20035,IFITM3,IFIT3,CXCL10,EPSTI1,SERPING1,LOC26010,OAS2,RSAD2,RTP4,"
test_genes1 = test_genes1.split(",")
test_genes2 = "NCKAP1L,THEMIS,HLA-DPB1,JUN,NFKBIA,PTPN6,TLR5,HLA-DRA,CD4,PRKCQ,PLCG1,CD247,TRBC1,FCGR1B,CCR7,TLR6,TLR10,FYB,TLR2,FCER1G,TRAT1,LILRB2,LILRB1,CLEC7A,BTLA,FYN,FGR,CD86,GRAP2,HLA-DQB1,ZAP70,PAG1,MNDA,IRAK3,LCK,DUSP6,PTPN22,PAK1,CTSH,CTLA4,LCP2,TLR4,CSK,FOS,LAT2,CARD11,HLA-DRB1,CHUK,CD24,UBASH3A,MAP2K6,TLR7,LY96,HCK,VASP,LYN,LIME1,MAPK14,MYD88,INPP5D,ITK,TRAC,SLA2,KLRK1,HLA-DMB,RPS6KA1,HLA-DQA1,MEF2C,SYK,MEF2A,TLR1,PTPRC,HLA-DOA,HLA-DOB,HLA-DMA,HLA-DPA1,CD19,CD3D,CD3E,BTK,CD3G"
test_genes2 = test_genes2.split(",")

# the species of interest:
species = "human"
# the GO terms aspect of interest:
aspect = ["biological_process"]
# ontologies to use:
target_onto = ["Gene Ontology", "Reactome", "Human Phenotype Ontology"]
# the genes of interest:
genes = test_genes2
# genes are symbols instead of UniProtKB IDS:
symbol = True

## data relative path:
data_path = "./data"

## URLs to download data files from and their respective resulting files:
# gaf file is downloaded automatically by a function call.
gaf_file = "%s/goa_%s.gaf" % (data_path, species)

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