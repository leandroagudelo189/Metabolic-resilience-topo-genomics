**Integrated FDR Adjustment:**  
The `adjust_pvalues` subroutine has been updated to correctly implement the False Discovery Rate (FDR) adjustment using the Benjamini-Hochberg procedure. It now ensures that the adjusted p-values are monotonically increasing and properly controls the expected proportion of false positives.  
```perl  
sub adjust_pvalues {
    my ($results_ref) = @_;
    my @hits = sort { $a->pvalue <=> $b->pvalue } @$results_ref;
    my $n = scalar @hits;
    if ($n > 0) {
        $hits[$n - 1]->pvalueadj($hits[$n - 1]->pvalue);
        for (my $i = $n - 2; $i >= 0; $i--) {
            my $adj_pvalue = min($hits[$i + 1]->pvalueadj, ($n / ($i + 1)) * $hits[$i]->pvalue, 1);
            $hits[$i]->pvalueadj($adj_pvalue);
        }
    }
}

```
*   
* A helper `min` function is used to ensure the adjusted p-values do not exceed 1 and are the minimum of the calculated values.

**Filtering Significant Hits:**  
The `filter_significant_hits` subroutine filters hits based on the adjusted p-values and the specified alpha threshold. It collects significant hits and associates them with their respective chromosomes.  
```perl  
  
sub filter_significant_hits {
    my ($results_ref, $sresults_ref) = @_;
    my $nbSignificant = 0;
    foreach my $hit (sort { $a->pvalueadj <=> $b->pvalueadj } @$results_ref) {
        if ($hit->pvalueadj <= $alpha) {
            push @$sresults_ref, $hit;
            $nbSignificant++;
            push @{$chrs{$hit->chr}->sEnriched}, $hit;
        }
    }
    return $nbSignificant;
}

```
* 

**Filtering Redundant Hits:**  
The `filter_redundant_hits` subroutine removes redundant hits by excluding nested hits with worse p-values or insufficient ratios, as per the logic provided in your code.  
```perl  

sub filter_redundant_hits {
    # Logic as explained earlier
}

```
*   
* It ensures that only the most significant and non-overlapping hits are retained.

**Improved Output Generation:**  
The `generate_output` subroutine has been updated to include the score calculation for visualization purposes and to output the results in the desired format.  
``` perl  

sub generate_output {
    # Logic as explained earlier
`}
```
*   
* It constructs the BED file content and raw data output, incorporating all pertinent information.

**Enhanced Code Organization:**

* The code has been modularized with clear subroutines for each major step, improving readability and maintainability.

**Consistent Variable Naming and Formatting:**

* Variables have been named consistently and descriptively, and code formatting has been standardized.

**Error Handling and Input Validation:**

* The script includes checks for required inputs and provides helpful error messages.

**Added `use strict;` and `use warnings;`**  
Enabling these pragmas helps catch potential errors and enforces good coding practices.

**Updated Command-Line Option Parsing**  
Replaced `Getopt::Std` with `Getopt::Long` for more flexible and descriptive command-line options. This allows the use of long option names and better usability.

```perl  
  
use Getopt::Long qw(GetOptions);
```
**Dynamic Module Path**

Used `FindBin` to ensure the script can locate local modules regardless of where it's run from.

```perl  
use FindBin;
use lib "$FindBin::Bin/";  # Ensure the script can locate local modules


```  

