# Sandbox

Directory of several simple scripts used for bioinformatic analysis

## genome2exome
Selects from a VCF files the variants overlaping the regions defined in a BED file plus X and Y bps upstream and downstream of each region, respectivelly

In the command line:
<pre><code>python genome2exome.py input.vcf regions.bed upstream_region_size downstream_region_size > output.vcf
</code></pre>
