# Introduction
Here, we will see some examples of how a user can apply TFTenricher to find relevant biological annotations to genes that are associated to a group of transcription factors (TFs)

## Basics
First, we import the TFTenricher file, which requires that the package has been properly installed. We recomment that a virtual environment is used, but this is optional. To install, run the following bash code when in the main folder:

```console
pip install .
```

Now, we can properly import and use TFTenricher


```python
# import TFTenricher
from TFTenricher import TFTenricher
```

Now, we start the analysis with some common immune-related TFs. Note that TFTenricher is designed to use gene symbols exclusively. 


```python
TFs = ['STAT1', 'GATA3', 'RELA', 'NFKB1', 'IRF4', 'STAT3', 'MYB']
enr = TFTenricher(TFs)
```

Now, we have mapped the TFs down to target genes using a correlation lookup table. In this example, we see how many of the immune-related TFs could not be found in the correlation table (0.0%). We also see that these seven TFs covered less than half of a percent of all TFs that were in the correlation table. 

We can view the putative target genes in the 'target_genes' attribute 


```python
print(enr.target_genes[:5])
print('Number of target genes: ', len(enr.target_genes))
```


Next, we can map these target genes to annotated gene sets, such as GO, KEGG, Reactome, the GWAS catalogue, or any set provided by the user. The default is GO. As per default, TFTenricher also performs a Benjamini-Hochberg multiple testing correction, but this can be omitted, or the user can provide a custom function.

This calculation can take around 30 seconds.


```python
enr.downstream_enrich()
enr.enrichments
```



Next, we want to plot the results.


```python
f, ax = enr.plot()
```

As can be seen, TFTenricher increases the power, and we thus get many significant pathways. We can therefore select to only plot the top N genes


```python
f, ax = enr.plot(plot_Ntop=5)
```

We can also sort on P-values


```python
f, ax = enr.plot(plot_Ntop=5, sorton='p')
```

## Optional inputs


### Providing a different multiple testing correction function
If the user wants an other method for correcting multiple testing, it can be provided in the downstream enrich process


```python
# Assuming that there already is a enr object of the TFTenricher class

# Print the number of terms passing FDR previously
print('Number of significant terms, Benjamini-Hochberg: ', enr.enrichments.FDR.sum())

# Define a Bonferroni correction (Note that the Bonferroni correction 
# is present in the stat_utils.py module, but we redefine it here for verbosity)
def bonferroni(p, FDR=0.05):
    return p < (FDR/len(p))


# Now, we plug in our Bonferroni function instead
enr.downstream_enrich(multiple_testing_correction=bonferroni)
print('Number of significant terms, Bonferroni: ', enr.enrichments.FDR.sum())

```

### Analysing other databases
The TFTenricher function also has several other databases to compare against, which is set by the 'db' option in the downstream_enrich method. The parameter 'db' can have the values 'REACTOME', 'KEGG', 'GWAS', 'GO', or be a dictionary set by the user.


```python
enr.downstream_enrich(db='GWAS')
enr.enrichments
```

### Comparing to a custom database of ontologies
To allow for swiching to a user-defined database, the 'db' parameter to the 'downstream_enrich' method should be a dictionary, with terms as keys and each set being a list/np.array of gene SYMBOL IDs


```python
# We define a custom ontology dict
my_own_db = {}

# we add a dummy ontology
# TFTenricher does not per default consider ontologies with less than 10 genes, so we add 50
my_own_db['pseudo_ontology'] = ['GENENAME' + str(i) for i in range(50)] 

# We also add one that highly matches the inferred target genes, for control
my_own_db['supermatch'] = enr.target_genes[:51]

# Now we check the results
enr.downstream_enrich(db=my_own_db)
enr.enrichments

```

### Testing other TF-target mappings than the correlation-based mapping
The default option of TFTenricher is to map TFs to target genes using a coexpression analysis. However, we can also map using any function that takes the form 'mapping_fun(TFs, silent, top_n_genes)', where TFs are the input list of TFs, 'silent' is whether to print output (bool), and top_n_genes is how many genes to include (int, set to None in the default -> TFTenricher performs a monte-carlo simulation for statistics on whether to include a downstream gene.   

#### Here, we set plug in the  built-in TRRUST mapping
NOTE: as presented in Sup. material 3, the TRRUST mapping is prone to biases, and we therefore recommend that the correlation-based approach is used. 


```python
from TFTenricher import map2trgt_utils

enr = TFTenricher(TFs, mapmethod=map2trgt_utils.trrust_genes)
enr.downstream_enrich()
enr.enrichments
```

#### We can also define our own mapping functions
The input paramater TFs is a list of TF SYMBOL names, while top_n_genes and silent are, as of now, hard-coded input variables. Here, we design a dummy mapping function that always returns the target genes GATA3 and GATA4.


```python
def foo(TFs, top_n_genes=None, silent=None):
    return ['GATA3', 'GATA4']


enr = TFTenricher(TFs, mapmethod=foo)
enr.target_genes

```

# Reference
The reference for the TFTenricher can be found by running the 'cite' method:


```python
enr.cite()
```

Note that 'cite()' also prints the references to the third-party resources that have been used in the enr object

COPYRIGHT (C) Rasmus Magnusson, 2020

Contact: rasma774@gmail.com
