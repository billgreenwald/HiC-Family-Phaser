# HiC-Family-Phaser

This script allows for the incorporation of family data with an already partially phased vcf to:
 1) Phase unphased variants where at least one member of the trio is phased
 2) Identify and fix switch errors according to inheritance.  This will not fix switch errors in founders
 3) Identify genotyping errors if a family is larger than a trio.
 
This script assumes that the input vcf has the majority of the inputs phased correctly, and should be phased from a source that does not utilize family structure (ie Hi-C data).

## Command Line Arguments:

### Required:

#### -i
The input vcf file to be phased

#### -o
The output file name

#### -s
The name of the sample within the vcf to phase

#### -p
Pedigree information file.  This file is a tab delimited file with 3 columns:  Child sample name, Paternal sample name, Maternal sample name

Example:
```
Child1    Father1    Mother1
Child2    Father2    Mother2
Child3    Father2    Mother2
```

In this example, Child2 and Child3 are siblings

### Optional:

#### -l
This flag indicates that a samples genotypes should not be changed during phase fixing.  The argument should be a comma separated list.  For example, if we did not want to change Father2's or Mother2's phasing, we would do:
```
-l Father2,Mother2
```

Please read the following section on recommended pipeline to see use of this argument.

## Recommended Pipeline

Running this script once will phase an inidividual and their parents to agree on haplotypes where they are unphased, fix switch errors when they are present, and unphase any genotyping errors.  Running this script on a single trio will therefore phase the family fine.  However, if you have a larger family, it will be necessary to run the script multiple times, and lock individuals that have already been phased so that they do not get changed to being unphased.  For example, given the family:
```
M1-----F1
    |
    |
    |
|-------|
|       |
C1      C2
```

We would first phase one trio, ie M1, F1, and C1.  This would create a vcf where M1, F1, and C1 all agree on phase.  Next, we want to incorporate any places where C2 has a switch error that agrees with the phasing from M1,F1, and C1, so we will phase M1, F1, and C2.  Finally, we want to flag all places where C1 and C2 disagree on phase within the family, so we will rephase M1, F1, and C1, but lock M1 and F1 so that their genotypes still agree with C2 after the script is done.

In a command line, this would be:
```
python fix_switch_errors.py -i input.vcf -o input-c1.vcf -s C1 -p structure.txt 
python fix_switch_errors.py -i input-c1.vcf -o input-c1-c2.vcf -s C2 -p structure.txt 
python fix_switch_errors.py -i input-c1-c2.vcf -o output.vcf -s C1 -p structure.txt -l M1,F1
```
