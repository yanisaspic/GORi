# GORi : Gene Ontology & Reactome itemsets

## About__________________________

*GORi* is a Python pipeline used to explore the relationships between genes and ontologies terms. The Gene Ontology (GO), the Reactome and the Human Phenotype Ontology (HPO) are three ontologies handled by the tool. Using *GORi* and a set of genes, you can determine which biological processes (GO), pathways (Reactome) and symptoms (HPO) your genes are involved in. Finally, *GORi* also highlights correlations between these biological elements (e.g. based on your genes set, are genes involved in the T-cell receptor pathway also involved in the autoimmunity ?)

## Setup__________________________

In order to use the pipeline, you first need to open the **settings.py** script located inside the main folder. 

#### Genes set
Indicate your target genes (*line.21*). You can either write them down in a Python list on the script directly (cf. the variable *test_genes*) or load them from a csv file (cf. *test_genes_3*). If you use a csv, make sure that :
1. there is no header.
2. the genes are separated by line breaks.

The genes can be identified using their *UniProtKB id* or their *symbol*. If you use genes symbols, make sure that the variable **symbol** (*l.25*) is set to ***True***.

#### GO annotations
A GO annotations file (GAF) can be downloaded from the Gene Ontology database. In this case, set the variable **download** to ***True*** (*l.27*) and attribute the file name to the variable **gaf_name** (*l.43*). You can navigate the available GAF files on the Gene Ontology website (http://current.geneontology.org/annotations/index.html).

If you prefer to use a local GAF file, set the **download** variable to ***False*** and indicate the path leading to your GAF file (*l.46*).

**Be careful:** you may need to comment (#) the lines 44 or 46 depending on whether you use a local GAF file or a downloaded one, respectively.

## Pipeline__________________________

The 3 steps of the pipeline are required for the tool to function correctly. None can be skipped.

### 1. run **download_files.py**

The Gene and Human Phenotype Ontologies are downloaded. Using files from the Reactome database, a Reactome ontology is also created.

After downloading the required files, a list of terms is associated for each existing gene found in the files. Each list includes all the leaf terms found using the annotations files, and all their parent terms, found using the ontologies.

At the end of this step, the downloaded files and a new **genes_annotations.json** file are stored in the **data/** subfolder.

### 2. run **mine_itemsets.py**

Using the new file **genes_annotations.json**, your gene set, and an Information Content estimator *[Ontology-based information content computation, Sánchez et al. (2011), Knowledge Based Systems 24 (297-303)]*, each term's weight (how deep in the ontology and specific the term is) and frequency (how many genes of the set are associated to the term) are computed. The results are displayed in a plot.

![terms metrics](https://user-images.githubusercontent.com/81527286/135145263-5c096d6e-381b-4fd9-bb5e-935c6accc57b.png)

A threshold value for the weight and the frequency must be inputted in the command panel. All the terms with a weight and a frequency value superior to the threshold will be used for the following itemset mining. If you input a threshold value preceded with a minus sign '**-**', only the terms with values inferior to the threshold will be kept instead.

The frequent itemset mining calculation time can vary a lot based on your selection's size.

**Be careful:** while the Gene and Human Phenotype ontologies give satisfying results with the IC method used, it should be noted that the Reactome ontology is homemade for technical reasons. Therefore, the weight values of the Reactome terms are very different from the others.

After inputting your threshold values and filtering the terms of interest, an FPTree algorithm is used to determine how often each term is paired with another. The resulting **association rules** are used by the application, and two metrics are calculated for each rule : its **confidence** and its **lift**.

An **association rule** is made of two elements : its body (or antecedent) and its head (or consequent) :

![Rule](https://user-images.githubusercontent.com/81527286/135145795-b6a13132-9b88-4be5-8dda-82f9168fc257.png)

Here, X is the body and Y is the head. This rule can be read as "X's presence implies Y's presence". In order to determine how significant a rule is, metrics are used. The basic metric in itemset mining is the **support** of the frequent itemset:

![support](https://user-images.githubusercontent.com/81527286/135145640-6e7537fd-7f40-4fb5-9111-dea3073a89c5.png)

The support of X is defined as the proportion of transactions (i.e. genes) in the dataset which contains the itemset X (i.e. the ontology term). Its value varies between 0 (none of the gene is associated to the term X) and 1 (all the genes of the set are associated with the term X). Using the support of an itemset, its **confidence** can be calculated :

![confidence](https://user-images.githubusercontent.com/81527286/135145145-4f217c0b-e8b9-4bfd-a01e-676f3ac11938.png)

The confidence of the (X => Y) rule is the ratio of the support of X and Y together to the support of X alone. It varies between 0 (X is never found in an itemset with Y) and 1 (X is never found in an itemset without Y). Please note that the confidence is asymetrical, i.e. conf(X => Y) != conf(Y => X). Finally, using the confidence, the **lift** of the rule can be calculated.

![lift](https://user-images.githubusercontent.com/81527286/135145253-c8783a5f-ad07-4a20-818e-d5116cb409db.png)

The lift of the (X => Y) rule is the ratio of the observed support to the support of X and Y if they were independent. Its value indicates how X and Y are dependent of each other. Here are 3 examples to illustrate this :

Say our lift l=3 : l>1, which means that X and Y are paired together 3 times more than what would be expected randomly. Perhaps the genes implied in the biological process X also play a part in the symptom Y ?

Now l=0.5 : l<1, which means that X and Y are paired together 2 times less than what would be expected randomly. Perhaps the genes responsible for the symptom Y are inhibited by the activity of the genes responsible for the biological process X ?

Finally, if l=1, the X and Y pair could be random and therefore the rule is of little interest.

Please note that biological interpretations of this metric are proposed, but the conclusions drawn from this analysis obviously depend on multiple factors, such as your gene set, your hypothesis, etc. Also, the lift is symetrical, i.e. lift(X => Y) == lift(Y => X).

The confidence and the lift are calculated for each pair of terms found using the FPTree algorithm. The two terms of each pair are from different ontologies, and the resulting metrics are plotted again :

![rules metrics](https://user-images.githubusercontent.com/81527286/135145326-759e7154-f110-4513-900d-c6de2b954e7f.png)

Each point on the plot represents an association rule of two terms. Once again, two threshold values are required, for the confidence and the lift of the rules of interest respectively.

The display application loading times can vary a lot based on your selection size.

**Be careful:** the lift value varies according the confidence value, and the confidence is dependent of the support. The use of two support-dependent metrics is debatable.

At the end of this second step, two new files **items.csv** and **rules.csv** are created in the subfolder **/results**.

### 3. run **app.py**

Using the Dash package, a Flask app is created in order to explore the data using a Cytoscape network. Once the application is ready, a local link is displayed in the command panel. Following this link leads you to the GORi webpage :

![GORi](https://user-images.githubusercontent.com/81527286/135145181-5d06910d-eb08-4c56-98d6-be203c094ee7.png)

The application can be split into 3 elements : a main network, a closeup network and a control panel.

1. the **main network** (center)

On this network, each term is represented by a node, and the rules involving a term are represented by the edges connected to its node.

The term's human readable-label is directly indicated on the corresponding node. The nodes are also size and color-coded :

![nodes' size and color](https://user-images.githubusercontent.com/81527286/135145292-682b7ba1-2c4f-40ff-88a9-a06cdf1c26b5.png) 

**Gene Ontology**, **Human Phenotype Ontology** and **Reactome** terms are indicated in **cyan**, **yellow** and **magenta** respectively. The size of the node reflects its frequency in the gene set : the bigger the node is, the more frequent the term is.

The color of the nodes also indicate their information content : the more specific a term is (= the higher its calculated weight), the brighter its node is :

![nodes' colors](https://user-images.githubusercontent.com/81527286/135145269-65fe0fce-ff70-463c-a8e5-c9e95379eb2f.png)

The edges are also color-coded. 

A red edge indicates a lift value lower than or equal to 0.5 : the terms are found together two times less than what could be expected if they were independent.

A green edge indicates the opposite : the terms are found together two times more than the expectations. In this case, the lift value is superior or equal to 2.

Finally, a black edge is used if 0.5 < lift < 2. These "indefinite" edges are hidden by default.

![edges](https://user-images.githubusercontent.com/81527286/135146543-d14818f2-00dc-4d9c-b570-f1c4ad359324.png)

You can move the entire network or a specific node around using your mouse. You can also zoom in and zoom out. The same can be done on the closeup network.

2. a **closeup network** (left)

This second network is similar to the main one, except it displays only a specific node and its edges. The specific node is the last node you've clicked on either networks. Using the closeup, you can directly navigate across the nodes of the network and focus on specific nodes of interest.

3. a **control panel** (right)

Using dropdowns and scrollbars, you can customize the information displayed on the networks.

The first dropdown allows you to rearrange the nodes according to different algorithms (from top left to bottom right : Klay, Euler, Concentric and Random).

![layout dropdown](https://user-images.githubusercontent.com/81527286/135145769-89cf5c8d-ae4e-453c-a59e-1686b65dceb7.png)
![layouts](https://user-images.githubusercontent.com/81527286/135145249-adc48fd9-33dc-49db-9c81-acc0b92568d0.png)

**Be careful:** the Klay algorithm can't handle a large amount of nodes and edges. Before using it, it might be necessary to filter the data of the network using the other panel tools.

The second dropdown uses keywords to quickly filter the displayed nodes and edges. The ontologies and the edges of interest are selected using the menu, and the network is updated according to the selected keywords :

![keywords](https://user-images.githubusercontent.com/81527286/135145202-7e5ef65b-5444-4836-ac71-ce183d302f35.png)

It is also possible to update the networks by using the metrics values, instead of global keywords. This method uses ranges of values : using the scrollbars corresponding to your metrics of interest, you can select ranges of interest and filter out all the elements outside of these ranges.

![ranges filter](https://user-images.githubusercontent.com/81527286/135145096-c12502f0-1ffa-4bf0-bedb-a235caded6f0.png)

Using the control panel, you can also view directly the metrics values corresponding of the selected node :

![node metrics](https://user-images.githubusercontent.com/81527286/135145256-6af2463f-f4a3-4297-9908-b51c7b84388b.png)

The connectivity metric corresponds to all the rules found involving the selected term. The term's label and ID are also displayed, and clicking them leads you to the corresponding webpage.
