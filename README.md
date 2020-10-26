# Linkage

Reference files required for running linkage on TKI server are located in /home/san/sbeecroft/Linkage_references

Further information on the linkdatagen software package can be found here: 
http://bioinf.wehi.edu.au/software/linkdatagen/

However, the version of linkdatagen we use is older than that on the website, so don't download any perl scripts from there.

Tp access executables from the Merlin package, you can either download the merlin software package from :
http://csg.sph.umich.edu/abecasis/merlin/download/ and then compile locally, and add to your PATH. otherwise, you can add my pre-exisiting merlin executables to your PATH (/home/san/SBeecroft/merlin/) OR you can place the pre-compiled executables found in this repo within your main linkage working directory. These are Merlin, Minx, pedmerge, pedwipe, pedstats. 

# Rationale
We use linkdatagen to manipulate genetic variants from either NGS/MPS data, or SNP-arrays. For NGS data, linkdatagen extracts informative variants from the sequencing data to use as datapoints for linkage analysis with Merlin.

# Usage
# To download:

_cd_ into your desired directory

_git clone https://github.com/SarahBeecroft/Linkage.git_

_cd Linkage_

_chmod 777 *.sh *.pl merlin minx pedwipe pedstats pedmerge_

_cp /home/san/sbeecroft/Linkage_references/* . _

# Before using:
You will need to create/upload the input files for the linkage program to use. The sequencing data input will be the bam files provided by the sequencing facility. These need to be uploaded into the Linkage working directory via filezilla or WinSCP. 

Other input files are simply text files- the familyID.ped, familyID.ws and familyID.vcflist.txt files, where 'familyID' is a prefix of your choosing. Spaces, slashes and other special characters are not allowed in the prefix! Make sure to name the files exactly as outlined here (case sensitive), otherwise the program will not work.  

These files are crucial to get right! Any small errror may confound your results without it being obvious. It's *vital* to keep the sample names within these files in a consistent order throughout. E.g. Dad, Mum, Bro, Sis, Bro2 in order throughout all the files. 

# Creating the familyID.ped file
This is a representation of the pedigree. The values must be separated by tabs or spaces.

Although pedigrees can become quite complex, all the information that is necessary to reconstruct individual relationships in a pedigree file can be summarized in five items: a family identifier, an individual identifier, a link to each parent (if available) and finally an indicator of each individual's sex. We then add whether an individual has the disease or not (affection status). 

As an example of how family relationships are described, we will construct a pedigree file for a small pedigree with two siblings, their parents and maternal grand-parents.

For this simple pedigree, the five key items take the following values:

FAMILY     PERSON   FATHER   MOTHER   SEX   AFFECTION_STATUS
example    granpa   unknown  unknown    m   unaffected
example    granny   unknown  unknown    f   affected
example    father   unknown  unknown    m   unaffected
example    mother   granpa   granny     f   affected
example    sister   father   mother     f   affected
example    brother  father   mother     m   affected

These key values constitute the columns of the pedigree file. Because of restrictions in early genetic programs, text identifiers are usually replaced by unique numeric values. 

Male = 1, Female = 2, Unknown = 0. Unaffected = 1, Affected = 2, Unknown = 0. 

In this pedigree, granpa and granny are what we call 'founders'- they are the first generation in the pedigree and therefore are allowed to have 'unknown' parents (i.e. mum and dad values of 0). father is also a founder, because he married into the family. All 'non-founders' require their parents to be specified by using the sample number of the mum or dad (in this case, mother has granpa and granny specified as her parents- 1 and 2 respectively

After replacing each identifier with unique identifier and recoding sexes/affection status, this is what a basic space-delimited pedigree file would look like:

<contents of basic.ped>
1   D18-0001   0  0  1  1
1   D18-0002   0  0  2  2
1   D18-0003   0  0  1  1
1   D18-0004   D18-0001  D18-0002  2  2
1   D18-0005   D18-0003  D18-0004  2  2
1   D18-0006   D18-0003  D18-0004  1  2
<end of basic.ped>

A pedigree file can include multiple families. Each family can have a unique structure, independent of other families in the dataset. 

*It's important to note that we need 

However, we will not always have sequencing data for everyone in the family! We need the .ws file to tell linkdatagen which samples we have data on, and which we don't. 

# Creating the familyID.ws file
The .ws file is a single line of N values, where N equals the number of individuals in the pedfile. The values must be separated by tabs or spaces. The kth value indicates where to find genotypes for the individual listed in the kth row of the pedigree file. The value  0 indicates that the corresponding individual is not genotyped. 

Using the above example, if we had sequencing data for granny, mother and sister only, the ws would be:
0 1 0 2 3 0

# Creating the familyID.vcflist.txt file
The vcflist.txt file is the full path of the .HM.vcf files that will be generated by linkdatagen. Basically, the format is 

/path/to/Linkage/working/dir/sampleID1.HM.vcf
/path/to/Linkage/working/dir/sampleID2.HM.vcf
/path/to/Linkage/working/dir/sampleID3.HM.vcf

Make sure the listed files are in the same order as the ped file and the ws file (obviously skipping any sampleIDs that didn't get sequenced). 

# Running samtools2merlin.sh
Usage is ./samtools2merlin.sh
The program will then ask for your input for:
1. A space seperated list of your sampleIDs of the bam files you are using (only include the sampleIDs of the files you have and want to use). 
2. The familyID prefix you have used for your ped, ws and vcflist files
3. The pattern of inheritance in the family- this can be dominant or recessive. 

The program will then start running. This may take some time. When it is finished, you will have two new directories. The one of interest is familyID.HapMap2_merlin, which has a sub-directory called genome. Within genome, the files you want are the parametric.pdf and parametric.tbl. 

