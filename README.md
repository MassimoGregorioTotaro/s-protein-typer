# S-protein-typer

## Description

The **S-protein-typer** application has been developed to provide support to the SARS-Cov2 barcoding typisation pipeline.
To analyze sequencing data obtained by the S1-ROI primer panel, files for analysis with the ARTIC pipeline are in the here included [archive](https://github.com/MassimoGregorioTotaro/s-protein-typer/blob/main/V1000SG.zip)
The original repository can be found on [GitLab](https://gitlab.com/MassimoGregorioTotaro/s-protein-typer.git)

## Visuals

![picture1](https://gitlab.com/MassimoGregorioTotaro/s-protein-typer/uploads/6ca30780d257eda4a55a8f66f108715b/pic1.png)

## Installation

The tool requires a [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment set up with the provided [environment file](https://github.com/MassimoGregorioTotaro/s-protein-typer/blob/main/s-protein-typer.yaml).
It has been tested on MacOS and GNU/Linux systems;
Windows 10 and above users are able to use it with minimal effort using [WSL](https://docs.microsoft.com/en-gb/windows/wsl/).

## Usage

A typical set up consists of:

- the 's_protein_typer.py' file itself;
- a 'reference' folder containing the reference FASTA files and, optionally, the classifier PKL file;
- a 'data' folder containing the FASTA files to be analysed;

The program can then simply be launched with:
**`$ python s_protein_typer.py`**;
it will analyse all the sequences in the 'data' folder with the default settings.

More options can be viewed running the program with the --help flag:
**`$ python s_protein_typer.py -h`**

- **--verbose** provides more information at runtime and for interpreting the output
- **--retrain** requires an aptly formatted file to retrain the provided classifier, the retrained classifier will then be exported in the same folder
- **--reference** is needed if the reference FASTA file is not in the default path
- **--slow** performs a Multiple Sequence Alignment, in case some sequences are particularly weak, it is more accurate, but also slower
- **--output** exports the alignment as a CSV file, to be quickly analysed again (e.g. with a different classifier) or used for retraining the classifier
- **--machine_readable** disables the terminal formatting so that the ouput cam be redirected or saved to a file via terminal
- **--classifier** in case the classifier is in a different path
- **--alignment** reads a pre-generated CSV file or an alignment FASTA file, instead of the sequences FASTA files, for a quicker classification analysis

## Support

Below you can find some commonly asked troubleshooting\usage tips.

The passed sequences are aligned both at the DNA level, to identify the correct ORF, and at the AA level, to identify the mutations. While the standard protocol, which involves a series of pairwise alignments between the reference and the target sequences, usually performs pretty well, in case the sequences are of poor quality or have extensive deletions around the start codon area, misalignments can occur. The user should normally be aware of a **messy output**, treat it with skepticism and thoroughly check the corrisponding sequence; one can reasonably suspect alignment issues in case of several reported INDEL mutations at the beginning of the sequence. To overcome the problem the calculation should be run again with the --slow flag (on subsets of the sequences dataset, if time is a strict concern, which must however comprise at least 5 'good' sequences per each problematic one).

The classification is performed by a pre-trained Random Forest classifier, provided as PKL file in the 'reference' folder.
Said **classifier** can, hovever, be **created anew and trained** by providing a correctly formatted CSV file (the original is provided [here](https://gitlab.com/MassimoGregorioTotaro/s-protein-typer/-/blob/main/reference/training_dataset.csv)) to the --retrain flag, which will also take care of saving it for future usage (beware, repeated runs will overwrite previous ones' outputs).
Being the generation and training steps non-deterministic, you might have to perform several runs with the -t flag until you get an optimal outcome.
The training dataset can be adjusted by adding/removing manually classified sequences: the mutations can be derived from a run of the s-protein-typer CSV output database, the class must be denoted in the first row, preceeding the sequence name with an underscore (e.g. B.1.1.7_XXX). By altering the training database one must be careful to provide good quality data, i.e. avoid incomplete sequences (the classifier will be underfit), do not provide too many duplicates (the classifier will be overfit), be consistent with class names (e.g. an identical sequence identified as both B.1.617.2_XXX and delta_XXX will cause serious classification issues), do not overlook curating the outgroup ('NA', it is fundamental for classfication accuracy) etc.

## Roadmap

Minor adjustments can be made, especially regarding the classifier, according to the development and evolution of the typisation efforts.

## Authors and acknowledgment

**Massimo G. Totaro**, Institute of Biochemistry, Graz University of Technology, Graz, Austria.

## License

[BSD-3-Clause](https://gitlab.com/MassimoGregorioTotaro/s-protein-typer/-/blob/main/)

## Project status

Maintained until end 2021.
