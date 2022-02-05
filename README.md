## Analyzed codes for the immune response of severe COVID-19 using scRNA-seq data

Here we submitted codes and methods on the manuscript entitled "Integrating single cell sequencing data with GWAS summary statistics reveals CD16+ monocytes and memory CD8+T cells involved in severe COVID-19".

# Abstract
Background: Understanding the host genetic architecture and viral immunity contributes to the development of effective vaccines and therapeutics for controlling the COVID-19 pandemic. Alterations of immune responses in peripheral blood mononuclear cells play a crucial role in the detrimental progression of COVID-19. However, the effects of host genetic factors on immune responses for severe COVID-19 remain largely unknown. 
Methods: We constructed a computational framework to characterize the host genetics that influence immune cell subpopulations for severe COVID-19 by integrating GWAS summary statistics (N = 969,689 samples) with four independent scRNA-seq datasets containing healthy controls and patients with mild, moderate, and severe symptom (N = 606,534 cells). We collected 10 pre-defined gene sets including inflammatory and cytokine genes to calculate cell state score for evaluating the immunological features of individual immune cells.
Results: We found that 34 risk genes were significantly associated with severe COVID-19, and the number of highly-expressed genes increased with the severity of COVID-19. Three cell-subtypes that are CD16+monocytes, megakaryocytes, and memory CD8+T cells were significantly enriched by COVID-19-related genetic association signals. Notably, three causal risk genes of CCR1, CXCR6, and ABO were highly expressed in these three cell types, respectively. CCR1+CD16+monocytes and ABO+ megakaryocytes with significantly up-regulated genes, including S100A12, S100A8, S100A9, and IFITM1, confer higher risk to the dysregulated immune response among severe patients. CXCR6+ memory CD8+ T cells exhibit a notable polyfunctionality including elevation of proliferation, migration, and chemotaxis. Moreover, we observed an increase in cell-cell interactions of both CCR1+ CD16+monocytes and CXCR6+ memory CD8+T cells in severe patients compared to normal controls among both PBMCs and lung tissues. The enhanced interactions of CXCR6+ memory CD8+T cells with epithelial cells facilitates the recruitment of this specific population of T cells to airways, promoting CD8+T cells mediated immunity against COVID-19 infection.
Conclusions: We uncover a major genetics-modulated immunological shift between mild and severe infection, including an elevated expression of genetics-risk genes, increase in inflammatory cytokines, and of functional immune cell subsets aggravating disease severity, which provides novel insights into parsing the host genetic determinants that influence peripheral immune cells in severe COVID-19.


# Introduction
Accumulating evidence have suggested alterations of immune responses in peripheral blood mononuclear cells (PBMCs) play a crucial role in the detrimental progression of COVID-19. A growing number of GWASs have identified numerous significant genetic variants associated with COVID-19 susceptibility and severity. Many earlier GWASs have shown that complex genetic dysregulations of peripheral immune cells with highly selective effects on the risk of immune-related diseases at the subcellular level. However, the effect of these genetic determinants on the peripheral immune cells for severe COVID-19 remains largely unknown. In light of there is no comprehensive study for revealing the genetically regulatory effects of peripheral immune cells on severe COVID-19, the present study is the first integrative genomic analysis by combining genetic information from GWAS with scRNA-seq data to genetically pinpoint immune cell types implicated in the etiology of severe COVID-19. 

# Methods
## 1. Single cell RNA-seq data on severe COVID-19  
In the current study, we downloaded four independent scRNA-seq datasets on COVID-19 and its severity in PBMC and BALF from the ArrayExpress database (Dataset #1: the accession number is E-MTAB-9357) from Su et al. study [11], and the Gene Expression Omnibus (GEO) database (Dataset #2: the accession number is GSE149689 from Lee et al. study [20], Dataset #3: the accession number is GSE150861 from Guo et al. study [12], and Dataset #4: the accession number is GSE158055 [9]). For dataset #1, this dataset contained 270 peripheral blood samples including 254 samples with different COVID-19 severity (i.e., mild N = 109, moderate N = 102, and severe N = 50) and 16 healthy controls for scRNA-seq analysis. For the dataset #2, there were eight patients with COVID-19 of varying clinical severity, including asymptomatic, mild, and severe, and four healthy controls with PBMCs. As for the dataset #3, there were five peripheral blood samples from two severe COVID-19 patients at three different time points during tocilizumab treatment, containing two different stages: severe stage and remission stage. With regard to the dataset #4, there were 12 BALF samples including three moderate and nine severe patients collected from lung tissues. For all datasets, the sample collection process underwent Institutional Review Board review and approval at the institutions where samples were originally collected. The COVID-19 severity was qualified by using the World Health Organization (WHO) ordinal scale (WOS), the National Early Warning Score (NEWS), or the Diagnosis and Treatment of COVID-19 (Trail Version 6). Single-cell transcriptomes for these four datasets were gathered by using the 10× Genomics scRNA-seq platform. 

## 2. GWAS summary statistics from the COVID-19 Host Genetic Consortium
  The meta-GWAS summary data on severe COVID-19 round 4 (B2_ALL, Susceptibility [Hospitalized COVID-19 vs. Population]) were downloaded from the official website of the COVID-19 Host Genetic Consortium [23] (https://www.covid19hg.org/; analyzed file named: “COVID19_HGI_B2_ALL_leave_23andme_20201020.txt.gz”; released date of October 4 2020). There were 7,885 hospitalized COVID-19 patients and 961,804 control participants from 21 independent contributing studies. The vast majority of participants in these contributing studies were of European ancestry (93%). The meta-GWAS summary statistics contained P values, Wald statistic, inverse-variance meta-analyzed log Odds Ratio (OR) and related standard errors. The 1,000 Genomes Project European Phase 3 [37] were used as a panel for pruning. Results from 23&Me cohort GWAS summary statistics were excluded from our current analysis. By filtering genetic variants without RefSNP number in the Human Genome reference builds 37, there were 9,368,170 genetic variants included with a major allele frequency (MAF) threshold of 0.0001 and the imputation score filter of 0.6. We used the qqman R package for figuring the Manhattan plot to visualize the meta-GWAS analysis results. The web-based software of LocusZoom [38] was utilized to visualize the regional association plots for identified risk loci (http://locuszoom.sph.umich.edu/).
  
 ## 3. Scripts:
  In the present sutyd, we leveraged numerous bioinformatics tools: linux-based tools incluidng MAGMA, S-MultiXcan, R-based tools including Rolypoly and permutation, and web-access tools inclduing the WEB-based Gene SeT AnaLysis Toolkit (WebGestalt; http://www.webgestalt.org) [42], the PhenoScanner V2 (http://www.phenoscanner.medschl.cam.ac.uk/) [45],  the Open Target Genetics (OTG, https://genetics.opentargets.org/) [46], STRING(v11.0, https://string-db.org/)[51], STITCH (v5.0, http://stitch.embl.de/)[53],ChEMBL (v2.6, https://www.ebi.ac.uk/chembl/) [54], and DGIdb database (https://www.dgidb.org/druggable_gene_categories). 
  In order to ensure our peers could follow our analyses, we have deposited the codes and methods in the current github, as the following example:
```
#compute rolypoly
######################
library("rolypoly")
library("data.table")
#
index<-c("normal","mild","moderate","severe")
lapply(index,function(x){
  file_n<-paste0("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_",x,"_cell.txt")
  merge_scexpr<-read.delim(file_n,sep = " ")
  colnames(merge_scexpr)<-annotation$V2
  merge_scexpr<-merge_scexpr[apply(merge_scexpr,1,sum)!=0,]
  #create the annotation files
  gene_name<-intersect(rownames(merge_scexpr),geneid_df1$label)
  geneid_df1<-geneid_df1[geneid_df1$label %in% gene_name,]
  merge_scexpr<-merge_scexpr[gene_name,]
  geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
#############################################
  file_na<-paste0("roly_",x,"_pre.RData")
  save(geneid_df1,merge_scexpr,file=file_na)
  })

 ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"

#sim_block_annotation$label<-rownames(merge_scexprc2)[1:1000]
#Rploy_remission_GSE.txt
 rolypoly_result <- rolypoly_roll(
   gwas_data = COVID19_GWAS_autosomes_maf2,
   block_annotation = geneid_df1,
   block_data = merge_scexpr,
   ld_folder =ld_path,
   bootstrap_iters = 100
  )
  save(rolypoly_result,file = "/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/rolypoly_mild_cell.RData")
```



# Reference
1.	Dong E, Du H, Gardner L: An interactive web-based dashboard to track COVID-19 in real time. Lancet Infect Dis 2020, 20:533-534.
2.	Wu Z, McGoogan JM: Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China: Summary of a Report of 72 314 Cases From the Chinese Center for Disease Control and Prevention. JAMA 2020.
3.	Berlin DA, Gulick RM, Martinez FJ: Severe Covid-19. N Engl J Med 2020.
4.	Richardson S, Hirsch JS, Narasimhan M, Crawford JM, McGinn T, Davidson KW, Barnaby DP, Becker LB, Chelico JD, Cohen SL, et al: Presenting Characteristics, Comorbidities, and Outcomes Among 5700 Patients Hospitalized With COVID-19 in the New York City Area. JAMA 2020, 323:2052-2059.
5.	Guan WJ, Ni ZY, Hu Y, Liang WH, Ou CQ, He JX, Liu L, Shan H, Lei CL, Hui DSC, et al: Clinical Characteristics of Coronavirus Disease 2019 in China. N Engl J Med 2020, 382:1708-1720.
6.	Xu L, Ma Y, Yuan J, Zhang Y, Wang H, Zhang G, Tu C, Lu X, Li J, Xiong Y, et al: COVID-19 Quarantine Reveals That Behavioral Changes Have an Effect on Myopia Progression. Ophthalmology 2021.
7.	Pedersen SF, Ho YC: SARS-CoV-2: a storm is raging. J Clin Invest 2020, 130:2202-2205.
8.	Takahashi T, Ellingson MK, Wong P, Israelow B, Lucas C, Klein J, Silva J, Mao T, Oh JE, Tokuyama M, et al: Sex differences in immune responses that underlie COVID-19 disease outcomes. Nature 2020, 588:315-320.
9.	Chen G, Wu D, Guo W, Cao Y, Huang D, Wang H, Wang T, Zhang X, Chen H, Yu H, et al: Clinical and immunological features of severe and moderate coronavirus disease 2019. J Clin Invest 2020, 130:2620-2629.
10.	Su Y, Chen D, Yuan D, Lausted C, Choi J, Dai CL, Voillet V, Duvvuri VR, Scherler K, Troisch P, et al: Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell 2020, 183:1479-1495.e1420.
11.	Guo C, Li B, Ma H, Wang X, Cai P, Yu Q, Zhu L, Jin L, Jiang C, Fang J, et al: Single-cell analysis of two severe COVID-19 patients reveals a monocyte-associated and tocilizumab-responding cytokine storm. Nat Commun 2020, 11:3924.
12.	Ren X, Wen W, Fan X, Hou W, Su B, Cai P, Li J, Liu Y, Tang F, Zhang F, et al: COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas. Cell 2021.
13.	Wen W, Su W, Tang H, Le W, Zhang X, Zheng Y, Liu X, Xie L, Li J, Ye J, et al: Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing. Cell Discov 2020, 6:31.
14.	Zhang JY, Wang XM, Xing X, Xu Z, Zhang C, Song JW, Fan X, Xia P, Fu JL, Wang SY, et al: Single-cell landscape of immunological responses in patients with COVID-19. Nat Immunol 2020, 21:1107-1118.
15.	Chua RL, Lukassen S, Trump S, Hennig BP, Wendisch D, Pott F, Debnath O, Thürmann L, Kurth F, Völker MT, et al: COVID-19 severity correlates with airway epithelium-immune cell interactions identified by single-cell analysis. Nat Biotechnol 2020, 38:970-979.
16.	Silvin A, Chapuis N, Dunsmore G, Goubet AG, Dubuisson A, Derosa L, Almire C, Hénon C, Kosmider O, Droin N, et al: Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19. Cell 2020, 182:1401-1418.e1418.
17.	Schulte-Schrepping J, Reusch N, Paclik D, Baßler K, Schlickeiser S, Zhang B, Krämer B, Krammer T, Brumhard S, Bonaguro L, et al: Severe COVID-19 Is Marked by a Dysregulated Myeloid Cell Compartment. Cell 2020, 182:1419-1440 e1423.
18.	Lee JS, Park S, Jeong HW, Ahn JY, Choi SJ, Lee H, Choi B, Nam SK, Sa M, Kwon JS, et al: Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19. Sci Immunol 2020, 5.
19.	Cao X: COVID-19: immunopathology and its implications for therapy. Nat Rev Immunol 2020, 20:269-270.
20.	Del Valle DM, Kim-Schulze S, Huang HH, Beckmann ND, Nirenberg S, Wang B, Lavin Y, Swartz TH, Madduri D, Stock A, et al: An inflammatory cytokine signature predicts COVID-19 severity and survival. Nat Med 2020, 26:1636-1643.
21.	Arunachalam PS, Wimmers F, Mok CKP, Perera R, Scott M, Hagan T, Sigal N, Feng Y, Bristow L, Tak-Yin Tsang O, et al: Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans. Science 2020, 369:1210-1220.
22.	The COVID-19 Host Genetics Initiative, a global initiative to elucidate the role of host genetic factors in susceptibility and severity of the SARS-CoV-2 virus pandemic. Eur J Hum Genet 2020, 28:715-718.
23.	Zhou S, Butler-Laporte G, Nakanishi T, Morrison DR, Afilalo J, Afilalo M, Laurent L, Pietzner M, Kerrison N, Zhao K, et al: A Neanderthal OAS1 isoform protects individuals of European ancestry against COVID-19 susceptibility and severity. Nat Med 2021, 27:659-667.
24.	Pairo-Castineira E, Clohisey S, Klaric L, Bretherick AD, Rawlik K, Pasko D, Walker S, Parkinson N, Fourman MH, Russell CD, et al: Genetic mechanisms of critical illness in COVID-19. Nature 2021, 591:92-98.
25.	Ma Y, Huang Y, Zhao S, Yao Y, Zhang Y, Qu J, Wu N, Su J: Integrative Genomics Analysis Reveals a 21q22.11 Locus Contributing Risk to COVID-19. Hum Mol Genet 2021.
26.	Gaziano L, Giambartolomei C, Pereira AC, Gaulton A, Posner DC, Swanson SA, Ho YL, Iyengar SK, Kosik NM, Vujkovic M, et al: Actionable druggable genome-wide Mendelian randomization identifies repurposing opportunities for COVID-19. Nat Med 2021, 27:668-676.
27.	Ellinghaus D, Degenhardt F, Bujanda L, Buti M, Albillos A, Invernizzi P, Fernández J, Prati D, Baselli G, Asselta R, et al: Genomewide Association Study of Severe Covid-19 with Respiratory Failure. N Engl J Med 2020, 383:1522-1534.
28.	Mapping the human genetic architecture of COVID-19. Nature 2021.
29.	Ren X, Wen W, Fan X, Hou W, Su B, Cai P, Li J, Liu Y, Tang F, Zhang F, et al: COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas. Cell 2021, 184:1895-1913.e1819.
30.	10x Genomics. https://www10xgenomicscom/solutions/single-cell/.
31.	Butler A, Hoffman P, Smibert P, Papalexi E, Satija R: Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol 2018, 36:411-420.
32.	Waltman L, van Eck NJ: A smart local moving algorithm for large-scale modularity-based community detection. The European Physical Journal B 2013, 86:471.
33.	Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh PR, Raychaudhuri S: Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 2019, 16:1289-1296.
34.	Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, Abecasis GR: A global reference for human genetic variation. Nature 2015, 526:68-74.
35.	Pruim RJ, Welch RP, Sanna S, Teslovich TM, Chines PS, Gliedt TP, Boehnke M, Abecasis GR, Willer CJ: LocusZoom: regional visualization of genome-wide association scan results. Bioinformatics 2010, 26:2336-2337.
36.	de Leeuw CA, Mooij JM, Heskes T, Posthuma D: MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol 2015, 11:e1004219.
37.	Wang J, Duncan D, Shi Z, Zhang B: WEB-based GEne SeT AnaLysis Toolkit (WebGestalt): update 2013. Nucleic Acids Res 2013, 41:W77-83.
38.	Kanehisa M, Goto S: KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res 2000, 28:27-30.
39.	Ma Y, Li J, Xu Y, Wang Y, Yao Y, Liu Q, Wang M, Zhao X, Fan R, Chen J, et al: Identification of 34 genes conferring genetic and pharmacological risk for the comorbidity of schizophrenia and smoking behaviors. Aging (Albany NY) 2020, 12:2169-2225.
40.	Ma Y, Qiu F, Deng C, Li J, Huang Y, Wu Z, Zhou Y, Zhang Y, Xiong Y, Yao J, Zhong Y, Qu J, Su J.: Analyzed codes for the immune response of severe COVID-19 using scRNA-seq data. https://githubcom/mayunlong89/COVID19_scRNA 2022.
41.	Barbeira AN, Dickinson SP, Bonazzola R, Zheng J, Wheeler HE, Torres JM, Torstenson ES, Shah KP, Garcia T, Edwards TL, et al: Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun 2018, 9:1825.
42.	Barbeira AN, Bonazzola R, Gamazon ER, Liang Y, Park Y, Kim-Hellmuth S, Wang G, Jiang Z, Zhou D, Hormozdiari F, et al: Exploiting the GTEx resources to decipher the mechanisms at GWAS loci. Genome Biol 2021, 22:49.
43.	Barbeira AN, Pividori M, Zheng J, Wheeler HE, Nicolae DL, Im HK: Integrating predicted transcriptome from multiple tissues improves association detection. PLoS Genet 2019, 15:e1007889.
44.	Ma X, Wang P, Xu G, Yu F, Ma Y: Integrative genomics analysis of various omics data and networks identify risk genes and variants vulnerable to childhood-onset asthma. BMC Med Genomics 2020, 13:123.
45.	Xu M, Li J, Xiao Z, Lou J, Pan X, Ma Y: Integrative genomics analysis identifies promising SNPs and genes implicated in tuberculosis risk based on multiple omics datasets. Aging (Albany NY) 2020, 12:19173-19220.
46.	von Mering C, Huynen M, Jaeggi D, Schmidt S, Bork P, Snel B: STRING: a database of predicted functional associations between proteins. Nucleic Acids Res 2003, 31:258-261.
47.	Szklarczyk D, Santos A, von Mering C, Jensen LJ, Bork P, Kuhn M: STITCH 5: augmenting protein-chemical interaction networks with tissue and affinity data. Nucleic Acids Res 2016, 44:D380-384.
48.	Cotto KC, Wagner AH, Feng YY, Kiwala S, Coffman AC, Spies G, Wollam A, Spies NC, Griffith OL, Griffith M: DGIdb 3.0: a redesign and expansion of the drug-gene interaction database. Nucleic Acids Res 2018, 46:D1068-d1073.
49.	Calderon D, Bhaskar A, Knowles DA, Golan D, Raj T, Fu AQ, Pritchard JK: Inferring Relevant Cell Types for Complex Traits by Using Single-Cell Gene Expression. Am J Hum Genet 2017, 101:686-699.
50.	Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ: Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience 2015, 4:7.
51.	Puram SV, Tirosh I, Parikh AS, Patel AP, Yizhak K, Gillespie S, Rodman C, Luo CL, Mroz EA, Emerick KS, et al: Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 2017, 171:1611-1624.e1624.
52.	Fajgenbaum DC, June CH: Cytokine Storm. N Engl J Med 2020, 383:2255-2273.
53.	Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan CH, Myung P, Plikus MV, Nie Q: Inference and analysis of cell-cell communication using CellChat. Nat Commun 2021, 12:1088.
54.	Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B: scCODA is a Bayesian model for compositional single-cell data analysis. Nat Commun 2021, 12:6876.
55.	Thomson W, Jabbari S, Taylor A, Arlt W, Smith D: Simultaneous parameter estimation and variable selection via the logit-normal continuous analogue of the spike-and-slab prior. Journal of the Royal Society Interface 2019, 16:20180572.
56.	Aitchison J: The statistical analysis of compositional data. Journal of the Royal Statistical Society: Series B (Methodological) 1982, 44:139-160.
57.	Ghoussaini M, Mountjoy E, Carmona M, Peat G, Schmidt EM, Hercules A, Fumis L, Miranda A, Carvalho-Silva D, Buniello A, et al: Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics. Nucleic Acids Res 2021, 49:D1311-d1320.
58.	Battle A, Brown CD, Engelhardt BE, Montgomery SB: Genetic effects on gene expression across human tissues. Nature 2017, 550:204-213.
59.	Wang Q, Chen R, Cheng F, Wei Q, Ji Y, Yang H, Zhong X, Tao R, Wen Z, Sutcliffe JS, et al: A Bayesian framework that integrates multi-omics data and gene networks predicts risk genes from schizophrenia GWAS data. Nat Neurosci 2019, 22:691-699.
60.	Ma Y, Li MD: Establishment of a Strong Link Between Smoking and Cancer Pathogenesis through DNA Methylation Analysis. Sci Rep 2017, 7:1811.
61.	Auwul MR, Rahman MR, Gov E, Shahjaman M, Moni MA: Bioinformatics and machine learning approach identifies potential drug targets and pathways in COVID-19. Brief Bioinform 2021.
62.	More SA, Patil AS, Sakle NS, Mokale SN: Network analysis and molecular mapping for SARS-CoV-2 to reveal drug targets and repurposing of clinically developed drugs. Virology 2021, 555:10-18.
63.	Bryois J, Skene NG, Hansen TF, Kogelman LJA, Watson HJ, Liu Z, Brueggeman L, Breen G, Bulik CM, Arenas E, et al: Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease. Nat Genet 2020, 52:482-493.
64.	Cortal A, Martignetti L, Six E, Rausell A: Gene signature extraction and cell identity recognition at the single-cell level with Cell-ID. Nat Biotechnol 2021, 39:1095-1102.
65.	Manne BK, Denorme F, Middleton EA, Portier I, Rowley JW, Stubben C, Petrey AC, Tolley ND, Guo L, Cody M, et al: Platelet gene expression and function in patients with COVID-19. Blood 2020, 136:1317-1329.
66.	Shaath H, Vishnubalaji R, Elkord E, Alajez NM: Single-Cell Transcriptome Analysis Highlights a Role for Neutrophils and Inflammatory Macrophages in the Pathogenesis of Severe COVID-19. Cells 2020, 9.
67.	Rydyznski Moderbacher C, Ramirez SI, Dan JM, Grifoni A, Hastie KM, Weiskopf D, Belanger S, Abbott RK, Kim C, Choi J, et al: Antigen-Specific Adaptive Immunity to SARS-CoV-2 in Acute COVID-19 and Associations with Age and Disease Severity. Cell 2020, 183:996-1012.e1019.
68.	King KR, Aguirre AD, Ye YX, Sun Y, Roh JD, Ng RP, Jr., Kohler RH, Arlauckas SP, Iwamoto Y, Savol A, et al: IRF3 and type I interferons fuel a fatal response to myocardial infarction. Nat Med 2017, 23:1481-1487.
69.	Hadjadj J, Yatim N, Barnabei L, Corneau A, Boussier J, Smith N, Péré H, Charbit B, Bondet V, Chenevier-Gobeaux C, et al: Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients. Science 2020, 369:718-724.
70.	Lippi G, Plebani M, Henry BM: Thrombocytopenia is associated with severe coronavirus disease 2019 (COVID-19) infections: A meta-analysis. Clin Chim Acta 2020, 506:145-148.
71.	Ma C, Cheung AF, Chodon T, Koya RC, Wu Z, Ng C, Avramis E, Cochran AJ, Witte ON, Baltimore D, et al: Multifunctional T-cell analyses to study response and progression in adoptive cell transfer immunotherapy. Cancer Discov 2013, 3:418-429.
72.	Akondy RS, Fitch M, Edupuganti S, Yang S, Kissick HT, Li KW, Youngblood BA, Abdelsamed HA, McGuire DJ, Cohen KW, et al: Origin and differentiation of human memory CD8 T cells after vaccination. Nature 2017, 552:362-367.
73.	Andrade F, Fellows E, Jenne DE, Rosen A, Young CS: Granzyme H destroys the function of critical adenoviral proteins required for viral DNA replication and granzyme B inhibition. Embo j 2007, 26:2148-2157.
74.	Li Y, Hou G, Zhou H, Wang Y, Tun HM, Zhu A, Zhao J, Xiao F, Lin S, Liu D, et al: Multi-platform omics analysis reveals molecular signature for COVID-19 pathogenesis, prognosis and drug target discovery. Signal Transduct Target Ther 2021, 6:155.
75.	Bruchez A, Sha K, Johnson J, Chen L, Stefani C, McConnell H, Gaucherand L, Prins R, Matreyek KA, Hume AJ, et al: MHC class II transactivator CIITA induces cell resistance to Ebola virus and SARS-like coronaviruses. Science 2020, 370:241-247.
76.	Wein AN, McMaster SR, Takamura S, Dunbar PR, Cartwright EK, Hayward SL, McManus DT, Shimaoka T, Ueha S, Tsukui T, et al: CXCR6 regulates localization of tissue-resident memory CD8 T cells to the airways. J Exp Med 2019, 216:2748-2762.
77.	Takamura S, Kato S, Motozono C, Shimaoka T, Ueha S, Matsuo K, Miyauchi K, Masumoto T, Katsushima A, Nakayama T, et al: Interstitial-resident memory CD8(+) T cells sustain frontline epithelial memory in the lung. J Exp Med 2019, 216:2736-2747.
78.	Zhao J, Yang Y, Huang H, Li D, Gu D, Lu X, Zhang Z, Liu L, Liu T, Liu Y, et al: Relationship between the ABO Blood Group and the COVID-19 Susceptibility. Clin Infect Dis 2020.
79.	Klok FA, Kruip M, van der Meer NJM, Arbous MS, Gommers D, Kant KM, Kaptein FHJ, van Paassen J, Stals MAM, Huisman MV, Endeman H: Incidence of thrombotic complications in critically ill ICU patients with COVID-19. Thromb Res 2020, 191:145-147.
80.	Grillet F, Behr J, Calame P, Aubry S, Delabrousse E: Acute Pulmonary Embolism Associated with COVID-19 Pneumonia Detected with Pulmonary CT Angiography. Radiology 2020, 296:E186-e188.
81.	Poran A, Harjanto D, Malloy M, Arieta CM, Rothenberg DA, Lenkala D, van Buuren MM, Addona TA, Rooney MS, Srinivasan L, Gaynor RB: Sequence-based prediction of SARS-CoV-2 vaccine targets using a mass spectrometry-based bioinformatics predictor identifies immunogenic T cell epitopes. Genome Med 2020, 12:70.
82.	Soudja SM, Ruiz AL, Marie JC, Lauvau G: Inflammatory monocytes activate memory CD8(+) T and innate NK lymphocytes independent of cognate antigen during microbial pathogen invasion. Immunity 2012, 37:549-562.
83.	Ziegler-Heitbrock L: The CD14+ CD16+ blood monocytes: their role in infection and inflammation. J Leukoc Biol 2007, 81:584-592.
84.	Kawanaka N, Yamamura M, Aita T, Morita Y, Okamoto A, Kawashima M, Iwahashi M, Ueno A, Ohmoto Y, Makino H: CD14+,CD16+ blood monocytes and joint inflammation in rheumatoid arthritis. Arthritis Rheum 2002, 46:2578-2586.
85.	Hambleton S, Goodbourn S, Young DF, Dickinson P, Mohamad SM, Valappil M, McGovern N, Cant AJ, Hackett SJ, Ghazal P, et al: STAT2 deficiency and susceptibility to viral illness in humans. Proc Natl Acad Sci U S A 2013, 110:3053-3058.
86.	Samji T, Khanna KM: Understanding memory CD8(+) T cells. Immunol Lett 2017, 185:32-39.
87.	Nie Z, Hu G, Wei G, Cui K, Yamane A, Resch W, Wang R, Green DR, Tessarollo L, Casellas R, et al: c-Myc is a universal amplifier of expressed genes in lymphocytes and embryonic stem cells. Cell 2012, 151:68-79.
88.	Jouan Y, Guillon A, Gonzalez L, Perez Y, Boisseau C, Ehrmann S, Ferreira M, Daix T, Jeannet R, François B, et al: Phenotypical and functional alteration of unconventional T cells in severe COVID-19 patients. J Exp Med 2020, 217.
89.	Finucane HK, Reshef YA, Anttila V, Slowikowski K, Gusev A, Byrnes A, Gazal S, Loh PR, Lareau C, Shoresh N, et al: Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat Genet 2018, 50:621-629.
90.	Xiang B, Deng C, Qiu F, Li J, Li S, Zhang H, Lin X, Huang Y, Zhou Y, Su J, et al: Single cell sequencing analysis identifies genetics-modulated ORMDL3+ cholangiocytes having higher metabolic effects on primary biliary cholangitis. Journal of Nanobiotechnology 2021, 19:406.
91.	Lv Y, Huang Y, Xu X, Wang Z, Yu Y, Ma Y, Wu M: Integrated multi-omics data analysis identifies a novel genetics-risk gene of IRF4 associated with prognosis of oral cavity cancer. medRxiv; 2021.
92.	Dong Z, Ma Y, Zhou H, Shi L, Ye G, Yang L, Liu P, Zhou L: Integrated genomics analysis highlights important SNPs and genes implicated in moderate-to-severe asthma based on GWAS and eQTL datasets. BMC Pulm Med 2020, 20:270.
93.	Psychiatric genome-wide association study analyses implicate neuronal, immune and histone pathways. Nat Neurosci 2015, 18:199-209.
94.	Bulik-Sullivan BK, Loh PR, Finucane HK, Ripke S, Yang J, Patterson N, Daly MJ, Price AL, Neale BM: LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet 2015, 47:291-295.

