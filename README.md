# Linkage

Reference files required for running linkage on TKI server are located in /home/san/SBeecroft/Linkage_refs

Further information on the linkdatagen software package can be found here: 
http://bioinf.wehi.edu.au/software/linkdatagen/

However, the version of linkdatagen we use is older than that on the website, so don't download any perl scripts from there.

Tp access executables from the Merlin package, you can either download the merlin software package from :
http://csg.sph.umich.edu/abecasis/merlin/download/ and then compile locally, and add to your PATH. otherwise, you can add my pre-exisiting merlin executables to your PATH (/home/san/SBeecroft/merlin/). 

# Rationale
We use linkdatagen to manipulate genetic variants from either NGS/MPS data, or SNP-arrays. For NGS data, linkdatagen extracts informative variants from the sequencing data to use as datapoints for linkage analysis with Merlin. The first step is to run auto_samtools4linkdatagen.sh on your bam files from the sequencing provider. This step extracts the informative variants for each person, and puts them into a .HM.vcf (one for each person). The following step 
