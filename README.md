# The Metagenome-Assembled Genome Inventory for Children (MAGIC) facilitates early-life gut bacteriome and virome studies


## Abstract

Existing microbiota databases are biased towards adult samples, hampering accurate profiling of the infant gut microbiome. Here, we generated a **M**etagenome-**A**ssembled **G**enome **I**nventory for **C**hildren (**MAGIC**) from a large collection of bulk and viral-like particle-enriched metagenomes from 0-7 years of age, encompassing `3,299` prokaryotic and `139,624` viral species-level genomes, `8.5%` and `63.9%` of which are unique to MAGIC. MAGIC improves early-life microbiome profiling, with the greatest improvement in read mapping observed in Africans. We then identified `54` candidate keystone species, including several *Bifidobacterium spp.* and four phages, forming guilds that fluctuated in abundance with time. Their abundances were reduced in preterm infants and were associated with childhood allergies. By analyzing the *B. longum* pangenome, we found evidence of phage-mediated evolution and quorum sensing-related ecological adaptation. Together, the MAGIC database recovers genomes that enable characterization of dynamics of early-life microbiomes, identification of candidate keystone species, and strain-level study of target species.


## Files

Please download `MAGIC` dataset from here: https://zenodo.org/doi/10.5281/zenodo.10369093
        
        

| Filename                                           | Filesize    |  SHA256SUM                                                         |
| -------------------------------------------------- | ----------- | ------------------------------------------------------------------ |
| MAGIC_pMAGs.tar.gz                                 | `5.57 GB`   | `9a293ed062e038ba5f47d53dacdad067f09b557ff720745caf2c22766872464a` |
| MAGIC_vMAGs.tar.gz                                 | `18.66 GB`  | `d0a63a1c19defa268ee35b65e798125bbeb60af5c01a580c3262925d43ec0bff` |
| MAGIC_K2DB.tar.gz                                  | `21.53 GB`  | `db76b92a97bb55391a78199bf6de04ff9795a7722864d099fa05f9772f19e2a9` |
| DataS1-Annotations_of_MAGIC_pMAGs.xlsx             | `6.13 MB`   | `0de2ceefda96d546e958f4c5425ca92f832097b679ea3204e4634fede7078399` |
| DataS2-Clustering_of_pOTUs.xlsx                    | `531.79 KB` | `00d2ec972f52c8852878126a738ad7f6fba3b03b8759b520ef5e134106cd2897` |
| DataS3-Annotations_of_MAGIC_vMAGs.xlsx             | `63.03 MB`  | `8d7e686588175d951b4afa15a347ee99e40f0efb868d842e4862bde2a63a734b` |
| DataS4-Clustering_of_vOTUs.xlsx                    | `17.48 MB`  | `59254a35de0179fd4f29197222c29f6d267442b3f51545524bf4a49840f4c4f3` |
| DataS5-12-Annotations_of_MAGIC-Proteins.tsv.gz     | `1.17 GB`   | `e2221581bac3f09c7dba6af2273241792e0e87f7182dd0a3240c7e5609a102eb` |
| Proteins_in_pMAGs_vMAGs.tsv.gz                     | `840.86 MB` | `91625f705b134bf2f76e7ac11cdeb91d47236a8b32129c1d7c70d9a7b5273ceb` |
| MAGIC_proteins_95.fasta.gz                         | `2.01 GB`   | `b3803ff4ce59346bb0e52f0c1876b4e8fdee3005c2e15290ff0957819d87100a` |

After downloading, please verify the completeness of the data:

```bash
sha256sum -c SHA256SUM.txt
```

### MAGs genomic FASTA file

#### MAGs Folder structure

Each MAG was assigned a unique `9-digit` ID. The MAGs are stored within subfolders named with the first three, middle and last three digits of their IDs. For example, the sequence for a pMAG numbered as `000000001` is stored in `MAGIC_pMAGs/000/000/001/MAGIC_pMAG_000000001.fa`.

#### MAGIC_pMAGs

`MAGIC_pMAGs.tar.gz`: This is a compressed folder including fasta format files of total of `26352` strain-level prokaryotic (bacterial pr archaeal) metagenome-assembled genomes (pMAGs). After uncompressing, the folder structure will appear as follows:

```sh
$ tar -xzvf MAGIC_pMAGs.tar.gz

MAGIC_pMAGs/
MAGIC_pMAGs/000/
MAGIC_pMAGs/000/000/
MAGIC_pMAGs/000/000/001/
MAGIC_pMAGs/000/000/001/MAGIC_pMAG_000000001.fa
MAGIC_pMAGs/000/000/001/MAGIC_pMAG_000000001.fa.seqkit.stats.tsv
MAGIC_pMAGs/000/000/002/
MAGIC_pMAGs/000/000/002/MAGIC_pMAG_000000002.fa
MAGIC_pMAGs/000/000/002/MAGIC_pMAG_000000002.fa.seqkit.stats.tsv
MAGIC_pMAGs/000/000/003/
MAGIC_pMAGs/000/000/003/MAGIC_pMAG_000000003.fa
MAGIC_pMAGs/000/000/003/MAGIC_pMAG_000000003.fa.seqkit.stats.tsv
........
```

#### MAGIC_vMAGs

`MAGIC_vMAGs.tar.gz`: This is a compressed folder including fasta format files of total of `191646` strain level viral metagenome-assembled genomes (vMAGs). After uncompressing, the folder structure will appear as follows:

```sh
$ tar -xzvf MAGIC_vMAGs.tar.gz

MAGIC_vMAGs/
MAGIC_vMAGs/000/
MAGIC_vMAGs/000/000/
MAGIC_vMAGs/000/000/001/
MAGIC_vMAGs/000/000/001/MAGIC_vMAG_000000001.fa
MAGIC_vMAGs/000/000/001/MAGIC_vMAG_000000001.fa.seqkit.stats.tsv
MAGIC_vMAGs/000/000/002/
MAGIC_vMAGs/000/000/002/MAGIC_vMAG_000000002.fa
MAGIC_vMAGs/000/000/002/MAGIC_vMAG_000000002.fa.seqkit.stats.tsv
MAGIC_vMAGs/000/000/003/
MAGIC_vMAGs/000/000/003/MAGIC_vMAG_000000003.fa
MAGIC_vMAGs/000/000/003/MAGIC_vMAG_000000003.fa.seqkit.stats.tsv
........
```

### MAGs's annotations and clustering information

#### DataS1-Annotations_of_MAGIC_pMAGs.xlsx (`26,352 entries` * `31` columns)

| Field Name           | Description                                                                                                              |
|----------------------|--------------------------------------------------------------------------------------------------------------------------|
| pMAG_id              | ID of the pMAG (primary key)                                                                                             |
| pOTU_id              | ID of the pOTU                                                                                                           |
| pMAG                 | Original name of the pMAG                                                                                                |
| pOTU                 | Name of the representative species-level pMAG                                                                            |
| project_accession    | Project from which the MAG was generated                                                                                 |
| assembly_group       | Group in which contigs of the MAG were assembled. It is the same as the sample ID, or as the subject ID when multiple samples of the same subject were available |
| completeness         | Completeness estimated by CheckM                                                                                         |
| contamination        | Contamination estimated by CheckM                                                                                        |
| strain heterogeneity | Strain heterogeneity estimated by CheckM                                                                                 |
| MIMAG_quality_level  | Quality level of MAG based on the standards of the Minimum Information about a Metagenome-Assembled Genome (MIMAG)       |
| SGB_quality_level    | Quality level of MAG based on the criteria of species level genomic bins (SGB)                                           |
| quality_score        | Completeness - 5 * contamination                                                                                         |
| classification       | Taxonomic assignment by GTDBtk                                                                                           |
| taxonomy             | Refined taxonomic assignment, used in the MAGIC database                                                                 |
| Length               | Length (bp) of the pMAG                                                                                                  |
| Count                | Number of contigs for the pMAG                                                                                           |
| GC (%)               | GC content (%) of the pMAG                                                                                               |
| N50                  | Length of the shortest contig for which longer and equal length contigs cover at least 50 % of the assembly              |
| pOTU_unique          | Uniqueness of the pOTU compared to the publicly available human gut pOTUs (yes: unique; no: overlapped with known pOTUs) |
| GUNC-n_genes_called  | Number of genes called by Prodigal                                                                                       |
| GUNC-n_genes_mapped  | Number of genes mapped by diamond into GUNC refDB                                                                        |
| GUNC-n_contigs       | Number of contigs containing mapped genes                                                                                |
| GUNC-taxonomic_level | Taxonomic clade labels at this taxonomic level were used to calculate values in all following columns. For each genome, all scores at six levels (species level can be added using a command-line option) are calculated  |
| GUNC-proportion_genes_retained_in_major_clades | Only major clades that have >2% of all mapped genes assigned to them are retained to calculate other scores. Value of this column is n_genes_retained/n_genes_mapped    |
| GUNC-genes_retained_index | n_genes_mapped/n_genes_called * proportion_genes_retained_in_major_clades, i.e. a portion of all called genes retained in major clades                                                       |
| GUNC-clade_separation_score | A result of applying a formula explained in GUNC paper to taxonomy and contig labels of genes retained in major clades. Ranges from 0 to 1 and is set to 0 when genes_retained index is <0.4 because that is too few genes left |
| GUNC-contamination_portion | Portion of genes retained in major clades assigned to all clades except the one clade with the highest proportion of genes assigned to it                                                   |
| GUNC-n_effective_surplus_clades | Inverse Simpson Index of fractions of all clades - 1 (as 1 genome is expected) describing the extent of chimerism, i.e. the effective number of surplus clades represented at a tax level |
| GUNC-mean_hit_identity              | Mean identity with which genes in abundant lineages (>2%) hit genes in the reference                                                                                               |
| GUNC-reference_representation_score | genes_retained_index * mean_hit_identity. Estimates how well a genome is represented in the GUNC DB                                                                                |
| GUNC-pass.GUNC                      | Overall assessment by GUNC. A genome passes if clade_separation_score <= 0.45, a cutoff benchmarked using simulated genomes                                                        |

#### DataS2-Clustering_of_pOTUs.xlsx

| Field Name    | Description                                                      |
| ------------- | ---------------------------------------------------------------- |
| MAGIC_pOTU_id | ID of the pOTU in MAGIC                                          |
| Rep_DB        | Source of the representative pOTU                                |
| Rep_FA        | File name of the representative pOTU                             |
| MAGIC         | List of pMAGs of the pOTU from MAGIC                             |
| CGR2          | List of pMAGs of the pOTU from CRG2                              |
| ELGG          | List of pMAGs of the pOTU from ELGG                              |
| GTDB          | List of pMAGs of the pOTU from GTDB                              |
| Hadza         | List of pMAGs of the pOTU from the Hadza hunter-gatherer dataset |
| IMGG          | List of pMAGs of the pOTU from IMGG                              |
| JMAG          | List of pMAGs of the pOTU from JMAG                              |
| SPMP          | List of pMAGs of the pOTU from SPMP                              |
| UHGG          | List of pMAGs of the pOTU from UHGG                              |
| WIS           | List of pMAGs of the pOTU from WIS                               |

#### DataS3-Annotations_of_MAGIC_vMAGs.xlsx (`191,646` entries * `36` columns)

| Field Name              | Description                                                                                                             |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------|
| vMAG_id                 | ID of the vMAG (primary key)                                                                                            |
| vOTU_id                 | ID of the vOTU                                                                                                          |
| vMAG                    | Original name of the vMAG                                                                                               |
| vOTU                    | Name of the representative species-level vMAG                                                                           |
| project_accession       | Project from which the MAG was generated                                                                                |
| assembly_group          | Group in which contigs of the MAG were assembled. It is the same as the sample ID, or as the subject ID when multiple samples of the same subject were available |
| viruses_type            | Type of the virus inferred by geNomad                                                                                   |
| contig_length           | Length (bp) of the vMAG                                                                                                 |
| provirus                | Existence of provirus determined by geNomad                                                                             |
| proviral_length         | Length of the provirus (bp) determined by geNomad                                                                       |
| gene_count              | Number of all genes determined by geNomad                                                                               |
| viral_genes             | Number of host genes determined by geNomad                                                                              |
| host_genes              | Number of host genes determined by geNomad                                                                              |
| checkv_quality          | Quality level calculated by CheckV                                                                                      |
| miuvig_quality          | Quality level of vMAG based on the standards of the Minimum Information about an Uncultivated Virus Genome (MIUVIG)     |
| completeness            | Completeness estimated by Checkv                                                                                        |
| completeness_method     | Method for the calculation of completeness used by CheckV                                                               |
| contamination           | Contamination estimated by Checkv                                                                                       |
| GC (%)                  | GC content (%) of the pMAG                                                                                              |
| N50                     | Length of the shortest contig for which longer and equal length contigs cover at least 50 % of the assembly             |
| taxonomy                | Refined taxaonomic assignemnt by geNomad and clustering, used in the MAGIC database                                     |
| species_all             | All species-level host predicted by iPHoP and Virus-Host-DB                                                             |
| species_best            | The best species-level hosts predicted by iPHoP and Virus-Host-DB                                                       |
| species_lca             | The Lowest Common Ancestor corresponding to all predicted species-level hosts                                           |
| species_lca_level       | Taxonomic level of the species_lca                                                                                      |
| host_phylum_best        | The best phylum-level hosts predicted by iPHoP and Virus-Host-DB                                                        |
| host_genus_all          | All genus-level hosts predicted by iPHoP and Virus-Host-DB                                                              |
| host_genus_best         | The best genus-level host predicted by iPHoP and Virus-Host-DB                                                          |
| host_genus_lca          | The Lowest Common Ancestor corresponding to all predicted genus-level hosts                                             |
| host_genus_best_lineage | Refined lineage of the host_genus_best, used in the calculation of virus-microbe-ratio                                  |
| PhaTyp_prediction       | Life style predicted by PhaTYP                                                                                          |
| PhaTyp_score            | Score of the PhaTyp prediction                                                                                          |
| public_vOTU             | Clustering with publicly available vOTU(s)                                                                              |
| size_vOTU               | Number of vMAGs in the vOTU                                                                                             |
| size_MAGIC              | Number of MAGIC-derived vMAGs in the vOTU                                                                               |
| vOTU_unique             | Uniquness of the vOTU compared to the publicly available human gut vOTUs (yes: unique; no: overlapped with known vOTUs) |

#### DataS4-Clustering_of_vOTUs.xlsx

| Field Name    | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| MAGIC_vOTU_id | ID of the pOTU in MAGIC                                      |
| Rep_DB        | Source of the representative vOTU                            |
| Rep_FA        | File name of the representative vOTU                         |
| MAGIC         | List of vMAGs of the vOTU from MAGIC                         |
| ELGV          | List of vMAGs of the vOTU from ELGV                          |
| GPD           | List of vMAGs of the vOTU from GPD                           |
| MGV           | List of vMAGs of the vOTU from MGV                           |
| GVD           | List of vMAGs of the vOTU from GVD                           |
| IMG_VR        | List of vMAGs of the vOTU from IMG_VR                        |
| RefSeq        | List of vMAGs of the vOTU from RefSeq                        |
| COPSAC_V      | List of vMAGs of the vOTU from COPSAC_V                      |
| JVD           | List of vMAGs of the vOTU from JVD                           |
| JP4D          | List of vMAGs of the vOTU from JP4D                          |
| Centenarians  | List of vMAGs of the vOTU from Centenarians                  |
| Hadza         | List of vMAGs of the vOTU from Hadza hunter-gatherer dataset |
| LOU           | List of vMAGs of the vOTU from LOU                           |
| HEVC          | List of vMAGs of the vOTU from HEVC                          |
| LLNEXT        | List of vMAGs of the vOTU from LLNEXT                        |


### MAGs's proteins and annotations

This table has three parts. The first part lists proteins in the pMAGs and vMAGs (`Proteins_in_pMAGs_vMAGs.tsv.gz`), whereas the second part provides functional annotations of the non-redundant proteins (`Data_S4-12_Annotations_of_MAGIC-Proteins.tsv.gz `). The third part provides the FASTA format sequence of each protein. Users may retrieve the list of genes on a MAG of interest and subsequently refer to the gene annotation table for annotations. Conversely, users may retrieve a list of MAGs carrying the genes of interest.

#### Proteins_in_pMAGs_vMAGs.tsv.gz (`70,538,090` entries * `5` columns)

| Field Name        | Description                                                                                         |
| ------------------|---------------------------------------------------------------------------------------------------- |
| MAG_id            | ID of the pMAG/vMAG                                                                                 |
| OTU_id            | ID of the pOTU/vOTU                                                                                 |
| source_mag        | Original name of the pMAG/vMAG                                                                      |
| original_protein  | ID of the protein annotated in the MAG (primary key)                                                |
| pv_rep            | ID of the representative protein. This is the foreign key refering to the primary key of Table S4b  |

#### Data_S5-12_Annotations_of_MAGIC-Proteins.tsv.gz (`9,548,653` entries * `49` columns)

| Field Name               | Description                                                                                                                                 |
|--------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| ID                       | ID of the protein (primary key)                                                                                                             |
| eggNOG_eggNOG_OGs        | eggNOG orthologous group                                                                                                                    |
| eggNOG_COG_category      | Clusters of Orthologous Genes (COG) category                                                                                                |
| eggNOG_Description       | Description of the COG category                                                                                                             |
| eggNOG_Preferred_name    | Mapping of seed ortholog to gene name                                                                                                       |
| eggNOG_GOs               | Gene Ontologies (GO)                                                                                                                        |
| eggNOG_EC                | Enzyme Commission (EC) annotation                                                                                                           |
| eggNOG_KEGG_ko           | Kyoto Encyclopaedia of Genes and Genomes (KEGG) orthology (KO)                                                                              |
| eggNOG_KEGG_Pathway      | KEGG pathway                                                                                                                                |
| eggNOG_KEGG_Module       | KEGG module                                                                                                                                 |
| eggNOG_KEGG_Reaction     | KEGG reaction                                                                                                                               |
| eggNOG_KEGG_rclass       | Classification of the KEGG reaction                                                                                                         |
| eggNOG_BRITE             | KEGG BRIATE identifier (a collection of hierarchical classification systems capturing functional hierarchies of various biological objects) |
| eggNOG_KEGG_TC           | Transporter in the Transporter Classification Database                                                                                      |
| eggNOG_CAZy              | Carbohydrate-active enzymes (CAZymes) annotated by eggNOG mapper                                                                            |
| eggNOG_BiGG_Reaction     | Reaction in the BiGG knowledgebase                                                                                                          |
| eggNOG_PFAMs             | Protein in the Protein families database annotated by eggNOG mapper                                                                         |
| VOGs                     | Viral genes in the VOGDB                                                                                                                    |
| VOG_best                 | The best VOG                                                                                                                                |
| VOG_best_cat             | Category of the best VOG                                                                                                                    |
| VOG_best_anno            | Description of the best VOG                                                                                                                 |
| AcrDB                    | Best hit to the computationally predicted anti-CRISPR (Acr) and Acr-associated (Aca) operon database                                        |
| UniRef_ID                | ID of the best hit to the UniProt Reference Clusters                                                                                        |
| UniRef_anno              | Annotation of the best UniRef hit                                                                                                           |
| Pfam                     | Protein in the Protein families database annotated by hmmsearch                                                                             |
| Pfam_anno                | Annotation of the Pfam protein                                                                                                              |
| KOfam                    | KO family annotated by kofam_scan                                                                                                           |
| KO_anno                  | Annotation of the KOfam                                                                                                                     |
| CAZy                     | CAZymes annotated by diamond                                                                                                                |
| SARG_sseqid              | Sequence ID in the Structured Antibiotic Resistance Gene (SARG) database                                                                    |
| SARG_Tag                 | Tag of the SARG (e.g., mutation, overexpression, regulator, repressor, etc)                                                                 |
| SARG_Type                | Type of antibiotic to which the SARG confers resistance (e.g., aminoglycoside)                                                              |
| SARG_Subtype             | Subtype of the SARG [e.g., aminoglycoside__AAC(3)-Ia]                                                                                       |
| SARG_HMM.category        | Name of the gene used as HMM profile [e.g., AAC(3)]                                                                                         |
| SARG_Mechanism.group     | Group of the mechanism of resistance of the SARG (e.g., Enzymatic inactivation)                                                             |
| SARG_Mechanism.subgroup  | Subgroup of the mechanism of resistance of the SARG (e.g., Acetyltransferases)                                                              |
| SARG_Mechanism.subgroup2 | Detail of the subgroup of the mechanism of resistance of the SARG [e.g., AAC(3)]                                                            |
| BRG_ID                   | ID of the biocide resistance gene (BRG) in the antibacterial Biocide & Metal Resistance Genes (BacMet) Database                             |
| BRG_Gene_name            | Name of the BRG                                                                                                                             |
| BRG_Compound             | Compound to which the BRG confers resistance                                                                                                |
| VFG                      | Virulence Factor gene (VFG) in the Virulence Factor Database (VFDB)                                                                         |
| VF_Name                  | Short name of the VFG                                                                                                                       |
| VF_FullName              | Full name of the VFG                                                                                                                        |
| VFCID                    | Category ID of the virulence factor                                                                                                         |
| VFcategory               | Category of the virulence factor                                                                                                            |
| Ig_like_protein          | Hit to the highly immunogenic outer capsid (HOC) protein (Ig-like)                                                                          |
| uniq_shared              | Uniquness of the protein compared to the publicly available proteins (uniq: unique; shared: overlapped with known proteins)                 |
| source_stat              | Summary of the source of the protein, expressed as "P_count V_count". E.g., a protein found in one pMAG and two vMAGs is marked as "P1V2"   |
| source                   | Category of the source of the protein, either from pmag(s), vmag(s), or "both" (pmag and vmag)                                              |

#### MAGIC_proteins_95.fasta.gz

A non-redundant protein dataset, `7,573,080` proteins (`95%` amino acid identity, AAI) predicted from the pMAGs by Bakta (v1.8.1) and `4,588,724` proteins (`95%` AAI) predicted from the vMAGs by prodigal-gv were clustered into `9,548,653` proteins at `95%` AAI.


### Workflow

#### MAGIC databases used for taxonomic profiling

`MAGIC_K2DB.tar.gz`: This is a phanta-style Kraken2 databases used for microbiome profiling. After uncompressing, the folder structure will appear as follows:

```sh
$ tar -xzvf MAGIC_K2DB.tar.gz
MAGIC_K2DB/

# Kraken2 hash table and other information
MAGIC_K2DB/taxo.k2d
MAGIC_K2DB/opts.k2d
MAGIC_K2DB/hash.k2d
MAGIC_K2DB/seqid2taxid.map
MAGIC_K2DB/inspect.out

# Bracken needed
MAGIC_K2DB/database.kraken
MAGIC_K2DB/database100mers.kraken
MAGIC_K2DB/database100mers.kmer_distrib
MAGIC_K2DB/database150mers.kraken
MAGIC_K2DB/database150mers.kmer_distrib
MAGIC_K2DB/database200mers.kraken
MAGIC_K2DB/database200mers.kmer_distrib

# MAGs information
MAGIC_K2DB/library/
MAGIC_K2DB/library/species_genome_size.txt
MAGIC_K2DB/library/strain_genome_size.txt
MAGIC_K2DB/library/prelim_map.txt

# Taxonomy information
MAGIC_K2DB/taxonomy/
MAGIC_K2DB/taxonomy/taxid.map
MAGIC_K2DB/taxonomy/names.dmp
MAGIC_K2DB/taxonomy/nodes.dmp
MAGIC_K2DB/taxonomy/prelim_map.txt
```

#### A workflow for taxonomic profiling based on the MAGIC database

1. Prepare the Phanta workflow

```sh
$ git clone -b magic_db https://github.com/ohmeta/phanta
# in an activated environment
conda env update --file /full/path/to/git/clone/phanta/phanta_env.yaml
# Or create a new environment
conda env create --name phanta --file /full/path/to/git/clone/phanta/phanta_env.yaml
```

Then please follow the documentation on github to install other dependences software.

2. Create a folder for the project

```
mkdir -p profiling_test
cd profiling_test
```

3. Prepare a sample sheet file `samples.rmhost.tsv`. E.g.,

| #sample_id  |   fq1                                  | fq2                                    |
| ----------- | -------------------------------------- | -------------------------------------- |
| ERR525724   | /full/path/to/ERR525724.rmhost.1.fq.gz | /full/path/to/ERR525724.rmhost.2.fq.gz |
| ERR525732   | /full/path/to/ERR525732.rmhost.1.fq.gz | /full/path/to/ERR525732.rmhost.2.fq.gz |
| ERR525735   | /full/path/to/ERR525735.rmhost.1.fq.gz | /full/path/to/ERR525735.rmhost.2.fq.gz |

4. Update config.yaml like below:

```sh
$ cp /full/path/to/git/clone/phanta/config.yaml ./
$ vim config.yaml

# Specify paths and threshold in the config.yaml, such as:

pipeline_directory: /full/path/to/git/clone/phanta

# Sample file specifies sample names and names of files containing sample reads
# Format: Tab-delimited, three columns
# sample_name  read1_file  [read2_file]
# if paired end, all samples must be paired-end
# if single end, all samples must be single-end
# See example (samp_file.txt) in the testing folder
sample_file: /full/path/to/profiling_test/samples.rmhost.tsv

# In which directory should results be outputted?
outdir: /full/path/to/profiling_test/results

# please uncompress MAGIC_K2DB.tar.gz
database: /full/path/to/MAGIC_K2DB

# Specifications for step one - classification of metagenomic reads
confidence_threshold: 0.1 # increase to reduce false positives - range from 0-1
gzipped: True # True or False - are the read files gzipped?
class_mem_mb: 42768 # memory in MB - minimum is the size of the Kraken2 database - must be at least 32 GB for the default database
class_threads: 16 # see usage instructions - can increase if you have more threads available; no need to change if you have fewer
single_end_krak: False # change if you would like to use the integrated prophage detection postprocessing script

# Specifications for step two - filtering false positive species
# essentially - what fraction of a viral genome should be covered to consider it a true positive?
#cov_thresh_viral: 0.10
cov_thresh_viral: 0.20
# how many unique minimizers should be covered in a viral genome ""?
minimizer_thresh_viral: 0
# same for bacteria
#cov_thresh_bacterial: 0.01
cov_thresh_bacterial: 0.02
minimizer_thresh_bacterial: 0
# archaea, eukaryotes
#cov_thresh_arc: 0.01
cov_thresh_arc: 0.02
minimizer_thresh_arc: 0
#cov_thresh_euk: 0
cov_thresh_euk: 0.005
minimizer_thresh_euk: 0

# Speciications for step three - per-species abundance estimation
read_length: 100 # if you change this, make sure you have an appropriate Bracken database built for this read length
filter_thresh: 10 # do not assign reads to species X if < this number of reads were classified to it

# Delete intermediate files? Examples in testing/classification/intermediate
delete_intermediate: False # True or False
```

5. Run phanta workflow based on MAGIC database

```sh
snakemake \
    --snakefile /full/path/to/git/clone/phanta/Snakefile \
    --configfile ./config.yaml \
    --until all \
    --cores 128 \
    --jobs 8
```


## Contact

+ Hein M Tun (heintun@cuhk.edu.hk)
+ Ye Peng (yepeng@cuhk.edu.hk)
