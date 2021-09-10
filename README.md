# it-workflow
### Sfanos Laboratory | JHMI
This is a version-controlled repository intended to serve as a resources for supporting statistical code and processed microbiome data for Jones et. al. (2021, JOURNAL TITLE) "Incorporation of data from multiple hypervariable regions when analyzing bacterial 16S rRNA gene sequencing data."

* Top-level directory: `it-workflow/`
```
./it-workflow/
├── analysis                        <-- contains analysis directory output  
│   ├── P00-tx-validation                         
│   ├── P03-summarize-qc                      
│   ├── P05-clust-tree-cm-goods
│   ├── P06-filtered-convert-tx-class                        
│   ├── P07-filtered-sum-taxa-db4.0                     
│   ├── P08-filtered-sum-alpha
│   ├── P13-summarize-qc                      
│   ├── P15-clust-tree-cm-goods
│   ├── P16-convert-tx-class                        
│   ├── P17-sum-taxa-db4.0.pl                    
│   ├── P18-sum-alpha
│   ├── P20-heatmap-species                         
│   ├── P21-alpha-vis                       
│   ├── P22-pcoa
│   ├── P23-filtered-tax-class                      
│   ├── P31-alpha-vis                       
│   ├── P32-pcoa
│   ├── P33-filtered-tax-class 
│   ├── P34-diff-ab-fresh-frozen 
│   └── P35-correl-fresh-frozen
├── code                            <-- contains code for analysis
├── data                            <-- contains data that may be version controlled
│   ├── dbs                         <-- location of dada2 16S dbs
│   ├── final                       
│   ├── processed
│   └── raw
└── external                        <-- too large to version control

```
