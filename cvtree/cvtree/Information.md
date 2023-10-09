#  cvtree project
## What is the data
Each file is a `<bacteria_name>.faa` that for testing purposes, includes a list of genes in the bacteria.
* The `>gi` lines are meta data just describing the name of the gene and unique identifier
* The rest of the lines are the actual genes. 
    Each letter is one of 24 amino acids

## What we're doing
Our goal is to **determine which out of the 41 bacteria are more closely related**. 
* Bacteria are closely related if they have more similar genes
* Instead of comparing genes, we can use statistics on the make up of the genes - specifically on the **kmers that make up the genes** 

### Example
**File: AcMNPV.faa**
>gi|9627743|ref|NP_054030.1| protein tyrosine phosphatase [Autographa californica nucleopolyhedrovirus]
MFPARWHNYLQCGQVIKDSNLICFKTPLRPELFAYVTSEEDVWTAEQIVKQNPSIGAIIDLTNTSKYYDG
VHFLRAGLLYKKIQVPGQTLPPESIVQEFIDTVKEFTEKCPGMLVGVHCTHGINRTGYMVCRYLMHTLGI
APQEAIDRFEKARGHKIERQNYVQDLLI

* If we look at a kmer of **length 6**, we get **MFPAR** 
    It's like a word -- a subsequence of a much longer gene
* **A k-mer** is then a **word or subsequence of a gene** of **k length**

**The first eight 6-mers**
1. MFPARW
2.  FPARWH
3.   PARWHN
4.    ARWHNY
5.     RWHNYL
6.      WHNYLQ...

## tidy.cpp
** The main function we're parallelising is `CompareAllBacteria()`**
** Step 1: Profile this function**

### `CompareAllBacteria()`
* A high correlation means 2 bacteria are very similar
* A low correlation means 2 bacterai are very different



# Code Replacement
## fopen_s to fopen

### fopen_s

FILE* bacteria_file;
errno_t OK = fopen_s(&bacteria_file, filename, "r");

if (OK != 0)
{
    fprintf(stderr, "Error: failed to open file %s\n", filename);
    exit(1);
}

### fopen

FILE* bacteria_file;
bacteria_file = fopen(filename, "r");
if(bacteria_file == NULL)
{
    fprintf(stderr, "Error: failed to open file %s\n", filename);
    exit(1);
}


