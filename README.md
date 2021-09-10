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
│   ├── P00-tx-validation.sh    
│   ├── P01-demux.sh
│   ├── P02-dada2.sh
│   ├── P03-sum-denoising-stats.pl   
│   ├── P04-prep-sfanos-db4.0.sh
│   ├── P05-clust-tree-cm-goods.sh
│   ├── P06-convert-tx-class.sh  
│   ├── P06-filtered-convert-tx-class.sh       
│   ├── P07-filtered-sum-taxa-db4.0.pl   
│   ├── P07-sum-taxa-db4.0.pl    
│   ├── P08-filtered-sum-alpha
│   ├── P08-sum-alpha
│   ├── P13-sum-denoising-stats.pl                     
│   ├── P15-clust-tree-cm-goods.sh
│   ├── P16-convert-tx-class.sh                      
│   ├── P17-sum-taxa-db4.0.pl                    
│   ├── P18-sum-alpha.pl
│   ├── P20-heatmap-species.r                        
│   ├── P21-alpha-vis.R                      
│   ├── P22-pcoa.R
│   ├── P23-filtered-tax-class.R                    
│   ├── P31-alpha-vis.R                       
│   ├── P32-pcoa.R
│   ├── P33-filtered-tax-class.R 
│   ├── P34-diff-ab-fresh-frozen.r 
│   └── P35-correl-fresh-frozen.r
├── data                            <-- contains data that may be version controlled
│   ├── P05-clust-tree-cm-goods
│   ├── P15-clust-tree-cm-goods                     
│   ├── manifests
│   ├── atcc-metadata.tsv
│   ├── atcc-metadata.txt
│   ├── atcc-percent-abundances.xlsx
│   ├── filtered-percent-abundance.xlsx
│   ├── patient-metadata.tsv
│   ├── patient-metadata.txt
│   ├── patient-simple-ids.xlsx
│   └── sfanos.2020-06-19.xlsx
└── external                        <-- too large to version control

```
