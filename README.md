# GUESSmyLT
An efficient way to guess the library type of your RNA-Seq data.
_________

 <img align="center" src="https://github.com/NBISweden/GAAS/edit/master/annotation/CheatSheet/pictures/library_types.jpg"/>

# library prep methods:

| kit | Description | Paired | Stranded | Strand according to mRNA | Strand according to `first strand`|
| --- | --- | --- | --- | --- 
| TruSeq RNA Sample Prep kit  | | yes | No | | fr-unstranded |
| SMARTer ultralow RNA protocol | | yes | No | | fr-unstranded |
| All dUTP methods, NSR, NNSR | | yes | Yes | RF | fr-firststrand 
| TruSeq Stranded Total RNA Sample Prep Kit | | yes | Yes | RF | fr-firststrand 
| TruSeq Stranded mRNA Sample Prep Kit | | yes | Yes | RF | fr-firststrand 
| NEB Ultra Directional RNA Library Prep Kit | | yes | Yes | RF | fr-firststrand 
| Agilent SureSelect Strand-Specific | | yes | Yes | RF | fr-firststrand 
| Directional Illumina (Ligation) | | yes | Yes | FR | fr-secondstrand
| Standard SOLiD | | Yes | yes | FR | fr-secondstrand
| ScriptSeq v2 RNA-Seq Library Preparation Kit | | yes | Yes | FR | fr-secondstrand
| SMARTer Stranded Total RNA | | yes | Yes | FR | fr-secondstrand 
| Encore Complete RNA-Seq Library Systems  | | yes | Yes | FR | fr-secondstrand
| NuGEN SoLo  | | yes | Yes | FR | fr-secondstrand
| Illumina ScriptSeq| |  yes | Yes | FR | fr-secondstrand
| SOLiD mate-pair protocol | | | | | ff

--rf orientation are produced using the Illumina mate-pair protocol.?


# External ressources:

[https://chipster.csc.fi/manual/library-type-summary.html](https://chipster.csc.fi/manual/library-type-summary.html)
[https://galaxyproject.org/tutorials/rb_rnaseq/](https://galaxyproject.org/tutorials/rb_rnaseq/)
[http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html)
[https://sailfish.readthedocs.io/en/master/library_type.html](https://sailfish.readthedocs.io/en/master/library_type.html)
[https://rnaseq.uoregon.edu](https://rnaseq.uoregon.edu)
[https://www.researchgate.net/post/What_is_the_difference_between_strand-specific_and_not_strand-specific_RNA-seq_data](https://www.researchgate.net/post/What_is_the_difference_between_strand-specific_and_not_strand-specific_RNA-seq_data)
