## **Guide to Perform PGE Analysis**

### **1\. Ensure You Have the Necessary Environment**

Before running the script, make sure you have the following:

* **Perl Installed:** Ensure that Perl is installed on your system. You can verify this by running:

``` bash
perl -v
```
* **Required Perl Modules:** The script uses several Perl modules. Install them using `cpan` or `cpanm` if they aren't already installed.  
```  bash  
 
cpan Class::Struct Getopt::Long Pod::Usage FindBin Hypergeometric
```    
* *Note:* If `Hypergeometric` is a custom module, ensure it's available in your `lib` directory or adjust the `use lib` path accordingly.  
* **Reference Datasets:** Ensure that all required reference dataset files (e.g., `ensembl42.txt`, `Mouse430_2.na24.annot.csv`, etc.) are present in the specified `$dataDir` (`/data2/pge/data` by default). If they're located elsewhere, use the `-p` option to specify the correct path.

### **2\. Prepare Your Query File**

Create a text file containing the list of gene names you want to analyze. Each gene name should be on a separate line or separated by spaces.

**Example: `query_genes.txt`**

### **3\. Verify Gene Symbols and Identifiers**

Ensure that the gene symbols in your query file match those used in the reference datasets. Gene identifiers can vary (e.g., Ensembl IDs vs. gene symbols). If necessary, map your gene symbols to the appropriate identifiers used in your reference dataset.

### **4\. Run the Improved Perl Script**

Use the command-line interface to execute your Perl script with the necessary options.

**Basic Command:**

```bash  

perl pge2.0.pl -q query_genes.txt
```
**With Additional Options:** You can customize the analysis by providing additional command-line options as needed.

**Specify Reference Dataset:**  
```bash  
 
-r ensembl42
```
* 

**Set Alpha Threshold (e.g., 0.01 for 1%):**  
```bash  

-a 0.01
```
* 

**Set Redundancy Filter Ratio (e.g., 2.0):**  
```bash  
 
-s 2.0
```
* 

**Map to Gene Symbols:**  
```bash  
 
-m symbols
```
* 

**Restrict to a Specific Chromosome (e.g., chromosome 1):**  
```bash  
  
-c 1
```
* 

**Specify Custom Data Directory:**  
```bash  

`-p /path/to/custom/data`
```
* 

**Example Command with All Options:**

```bash  
 
perl pge_script.pl -q query_genes.txt -r ensembl42 -a 0.01 -s 2.0 -m symbols -c 1 -p /custom/data/directory
```
### **5\. Interpret the Output**

The script generates several sections in the output, encapsulated within XML-like tags. Here's what each section represents:

* **`<mapping>`:** Details the mapping between your query genes and the reference dataset.  
* **`<raw>`:** Contains raw data of significant enriched regions, including p-values and adjusted p-values.  
* **`<bed>`:** Provides output in BED format, which can be used for visualization in genome browsers like UCSC Genome Browser or IGV.  
* **`<alpha>`:** Displays the alpha threshold used for determining significance.  
* **`<adjust>`:** Indicates the method used for p-value adjustment (e.g., FDR).  
* **`<reference>`:** Specifies the reference dataset used for the analysis.

### **6\. Visualize the Enriched Regions**

To visualize the enriched regions identified by the script:

1. **Export BED File:**  
   * Save the `<bed>` section of the output to a file, e.g., `enriched_regions.bed`.  
2. **Use a Genome Browser:**  
   * **UCSC Genome Browser:**  
     * Navigate to the UCSC Genome Browser.  
     * Select the appropriate genome assembly (e.g., mm10 for mouse).  
     * Use the "Add Custom Tracks" feature to upload your `enriched_regions.bed` file.  
     * Visualize the enriched regions alongside other genomic annotations.  
   * **Integrative Genomics Viewer (IGV):**  
     * Download and install IGV.  
     * Load the mouse genome assembly.  
     * Import the `enriched_regions.bed` file to visualize the enriched regions.
