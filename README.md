# Workshop on mining of CRISPR-Cas systems in cold-adapted bacteria

**Prediction of Antiphage defense systems in cold-adapted bacteria, with focus on CRISPR-Cas systems**

Background:

Organisms that inhabit cold environments are commonly classified into two overlapping groups: psychrophiles and psychrotolerants (or psychrotrophs). Psychrophiles have an optimal growth temperature of around 15°C and a maximum growth temperature of 20°C, while psychrotolerants grow optimally around 20°C and have a maximum growth temperature of 30°C. Psychrophiles predominate in marine ecosystems, whereas bacteria isolated from cold terrestrial environments are most often found to be psychrotolerant. Here, we employ the term ‘Cold-adapted bacteria’ referring to both psychrophilic and psychrotolerant bacteria. 

Bacteria encode multiple lines of defense against phages. These include defense systems with well-studied mechanisms, such as abortive infection (Abi) systems, re-striction-modification (RM) systems and CRISPR-Cas systems, and many recently discovered systems with lesser-known or unknown modes of action. CRISPR-Cas system recognizes and targets foreign nucleic acids like RM but is distinguished by the ability to acquire new specificities. CRISPR-Cas confers defense against invading genetic elements by integrating short fragments of foreign DNA into the CRISPR array, termed spacers. This integration enables subsequent identification and degradation of complementary nucleic acids (protospacers).

CRISPR-Cas system has gained a lot of attention for genetic engineering and genome editing, owing to its programmable RNA-guided endonuclease activity. Classification of CRISPR-Cas systems is based on effector complexes, yielding two primary classes and six distinct types: Class 1 (including types I, III, and IV) and Class 2 (including types II, V, and VI). Furthermore, these types can be divided into at least 34 subtypes. A fundamental divergence between Class 1 and Class 2 lies in their utilization of multi-Cas effector complexes and single effector nucleases, respectively, which renders Class 2 systems particularly suitable for genome editing and genome engineering applications 

In this workshop, we map the prevalence and distribution of antiphase defence systems in cold-adapted bacteria from various cold environments

Materials: 

Cold-adapted bacterial genome sequences downloaded from five databases: BacDive (https://bacdive.dsmz.de/), TEMPURA (http://togodb.org/db/tempura), MarRef v.1.7 (https://sfb.mmp2.sigma2.no/), MarDB v.1.6 (https://sfb.mmp2.sigma2.no/) and Ocean Microbiomics Database (OMD) (https://microbiomics.io/ocean/) v.1.1 

Data Filtering:

* Genomes from bacteria isolated only above 60°N and below 50°S, where ocean surface temperatures <10°C
* Single amplified genomes (SAGs) were excluded 
* Metagenome-assembled genomes (MAGs) from MarDB and OMD databases, were quality assessed using CheckM 
* Low-quality sequences with low completeness (<90%) or high contamination (>5%), were excluded 
* All genomes were taxonomically reclassified using the GTDB-Tk (v2.3.0) workflow against the Genome Taxonomy Database (GTDB) r.214 and archaea removed

**938 genomes were downloaded from the five databases for downstream analysis**

Run Bioinformatics tool using following criteria:

* Run CRISPRCasTyper (cctyper v1.8.0) (https://crisprcastyper.crispr.dk/#/submit) at default parameters [subtype probability above 0.754]. Only CRISPRs part of an intact CRISPR-Cas loci include
* Predict other prokaryotic antiphage defense systems using Prokaryotic Antiphage Defense LOCator: PADLOC (v1.1.4) (https://padloc.otago.ac.nz/padloc/) with parameter E-value <0.01 and coverage >0.8
* Plot phylogenetic trees using phyloseq and ggtree R package. 
* Plot Principal component analysis (PCA) to examine possible correlations between the low prevalence of CRISPR-Cas systems against the high prevalence of other systems, such as dXTPases. PCA analysis is performed at the taxonomic rank "family" using the factoMineR v2.8 and factoextra v1.0.7 packages in R

CRISPR-Cas systems are less frequent in cold-adapted bacteria, compared to mesophilic and thermophilic species
