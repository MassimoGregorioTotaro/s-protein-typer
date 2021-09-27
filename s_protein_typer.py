from argparse                import ArgumentParser
from Bio.Align.Applications  import MuscleCommandline
from Bio.AlignIO             import read
from Bio.Seq                 import Seq
from Bio.SeqIO               import parse, write
from Bio.SeqRecord           import SeqRecord
from collections             import defaultdict
from dill                    import dump, load
from io                      import StringIO
from joblib                  import Parallel, delayed
from glob                    import glob
from pandas                  import concat, DataFrame, read_csv, Series
from re                      import match, search, split
from sklearn.ensemble        import RandomForestClassifier
from sklearn.model_selection import cross_validate
from typing                  import List, Optional, Tuple, Union

THRESHOLD, VERBOSE = None, False

    #   Terminal output formatting escape sequences and padding values
FORMATTING = {"ITA"   : '\033[3m'
             ,"BOL"   : '\033[1m'
             ,"G_FG"  : '\033[92m'
             ,"R_FG"  : '\033[91m'
             ,"Y_FG"  : '\033[93m'
             ,"OFF"   : '\033[0m'
             ,"SEQ_C" : '\033[3;30;42m'
             ,"VAR_C" : '\033[3;30;41m'
             ,"MUT_C" : '\033[3;30;43m'}

SEQ_L, VAR_L = 25, 11

    #   SARS-CoV2 S-protein sequence loading specific optimisation values
GEN_MIN, GEN_MAX, S_START, S_END, LEN_MIN = 28000, 30000, 21000, 25500, 1250

class Classifier:
    """The Classifier class provides all the functions needed to make predictions within the s_protein_typer workflow"""
        #   AA to int encoding for the classifier
    encoder = defaultdict(lambda : 1, zip(list("-=XRKEHDAGCVPLIMFWNQSTY"), range(-1,22)))


    def __init__(self, cls:RandomForestClassifier, ref:Seq, threshold:float=0.9):
        """Initialise the Classifier class; it requires a sklearn classifier instance, 
           the reference sequence and a threshold value for the confidence interval"""
        self.cls = cls
        self.ref = ref
            #   Create an index to only load the residue positions also present in the reference, i.e. no insertions
        self.indices = [str(x+1) for x in range(len(self.ref))]
        self.threshold = threshold


    def save(self, file_name:str) -> None:
        """Allow Classifier persistence by dumping a .pkl file"""
        with open(file_name, "wb") as f:
            dump(self, f)


    def readCSV(self, file_name:str) -> DataFrame:
        """Read the given .csv file and return it as a DataFrame"""
            #   The CSV is read and transposed excluding the reference sequence
        df = read_csv(file_name, index_col=0).T.iloc[1:]
        df.columns = df.columns.astype("str")
        return DataFrame(df, columns=self.indices).fillna('=')
        

    def train(self, df:DataFrame):
        """Train the classifier with the passed DataFrame using cross validation"""
            #   Classes are derived from the DataFrame indexes, like: B.1.1.7_000
        classes = [s.split('_')[0] for s in df.index]
        self.classes = list(set(classes))
            #   The X and y training data are numerically encoded
        X = df.replace(self.encoder)
        y = [self.classes.index(x) for x in classes]
            #   Cross validation is used to train the classifier algorithm
        cv = cross_validate(self.cls, X, y, cv=5, scoring="recall_macro", return_estimator=True)
            #   The best estimator is selected for predictions
        self.cls = cv["estimator"][cv["test_score"].argmax()]


    def predict(self, s:Series) -> str:
        """Return the class prediction for the given Series"""
        def eval_threshold(t) -> str:
            """Format the output according to the confidence interval values"""
            if (estimate := t > 0.5).any() and (out := self.classes[estimate.argmax()]) != "NA":
                return ((t >= self.threshold).any() and FORMATTING['BOL'] or FORMATTING['ITA']) + f"{out:<10}{FORMATTING['OFF']}"
            return "NA"
        encoded = [[self.encoder[s[i]] if i in s.index else 0 for i in self.indices]]
        return eval_threshold(self.cls.predict_proba(encoded)[0])


def silent_print(string: str) -> None:
    if VERBOSE:
        print(f"-- {string}")


def read_fasta(*files:str, optimise:bool=False) -> List[SeqRecord]:
    """Read all SeqRecords in the passed fasta files, while processing them."""
    out = list()
    for f in files:
        silent_print(f"Reading '{f}':")
        for sqn in parse(f, "fasta"):
            (seq := sqn).seq = sqn.seq.upper().strip("NX")
            if optimise:
                    #   In case of full genomic data, only extract the S-protein expected region
                if GEN_MIN < len(seq) < GEN_MAX:
                    seq = seq[S_START:S_END]
                    seq.id = '∆' + seq.id
                    #   If the sequence is too short, don't consider it at all, else the output will be messy
                elif len(seq) < LEN_MIN:
                    seq.seq = Seq('---')
                    seq.id = '-' + seq.id
            silent_print(f"\tsequence '{seq.id}' loaded;")
            out.append(seq)
    return out


def align_sequences(*record:SeqRecord, reference:Optional[SeqRecord]=None) -> List[SeqRecord]:
    """Align all the passed SeqRecords sequences."""
    write(record, (handle := StringIO()), "fasta")
    alignments = MuscleCommandline(diags=True)(stdin=handle.getvalue())[0]
    alignment = list(read(StringIO(alignments), "fasta"))
    if (ids := [s.id for s in alignment]):
            #   Empty sequences are not read from the alignment, so they must be reinserted manually
        for s in record:
            if s.id not in ids:
                alignment.append(s)
                alignment[-1].seq = Seq('-' * len(alignment[0].seq))
        #   Put the reference sequence at the top of the list
    if reference:
        ref_pos = [i for i,x in enumerate(alignment) if x.id == reference.id][0] or 0
        alignment[0], alignment[ref_pos] = alignment[ref_pos], alignment[0]
    return alignment


def translate_seq(seq:Seq) -> Seq:
    """Translate the given sequence returning the longest found AA chain in the current ORF."""
    return sorted((seq + 'NNN'[:-(len(seq) % 3)]).translate().split('*'), key=len, reverse=True)[0].strip('X')


def translate_DNA_alignment(DNA_alignment:List[SeqRecord]) -> None:
    """Translate in-place all Seqs in the SeqRecord List, considering the reference ORF."""
    start = DNA_alignment[0].seq.find("ATG")
    for rec in DNA_alignment:
        rec.seq = translate_seq(rec.seq[start:].ungap('-'))


def classify(classifier:Optional[Classifier], series:Series) -> str:
    """Given the Classifier, predict the variant from the given Series."""
    return classifier and classifier.predict(series) or "NA"


def mutations_string(ref:Series, ser:Series) -> str:
    """Print the difference string between the target and reference Series."""
    out, str_to_int = [""], (lambda s: int(search('[0-9]+', s).group()) if s else 0)
    for src,trg in zip(ref.items(), ser.items()):
            #   Indexes and aminoacids for reference and target
        (_, sv), (ti, tv) = src, trg
            #   Index and overextension of the previous mutation
        fst, snd = (split('-', out[-1]) * 2)[:2]
            #   Is the current mutation index immediately following the previous' one?
        adjacent = 0 <= str_to_int(ti) - str_to_int(snd) < 2
            #   Check whether there's a mutation and, if so, report it, making sure to
            #   group adjacent insertions, deletions or unrecognised residues together
        if tv != '=':
            if tv == 'X':
                if ('*' in fst) and adjacent:
                    out[-1] = f"{fst}-{ti}"
                else:
                    out.append(f"*{ti}")
            elif tv == '-':
                if sv == '-':
                    continue
                if ('∆' in fst) and adjacent:
                    out[-1] = f"{fst}-{ti}"
                else:
                    out.append(f"∆{ti}")
            elif sv == '-':
                if ('+' in fst) and adjacent:
                    out[-1] = f"{fst}{tv}"
                else:
                    out.append(f"+{str(ti).split('_')[0]}_{tv}")
            else:
                out.append(f"{sv}{ti}{tv}")
        #   Enforce the delta displaying convention
    out = delta_convention_enforcer(out)
    return ' '.join((f"{s:<8}" for s in out if s))


def get_differences(AA_alignment:Union[List[SeqRecord], DataFrame], classifier:Optional[Classifier]) -> Tuple[str, DataFrame]:
    """Return the difference string and DataFrame for the sequences in the given alignment."""
    def compare_AAs(AA_alignment:Union[List[SeqRecord], DataFrame]) -> DataFrame:
        """Codify the mutations in the alignment sequences into a DataFrame."""
        if type(AA_alignment) is DataFrame:
            return AA_alignment
        i_offset, i_length = 0, 0
        df_index, sequences = zip(*((rec.id, rec.seq) for rec in AA_alignment))
        df = DataFrame(index=df_index)
            #   Run through the aligned sequences in parallel, calculate the correct index and find the AA differences
        for i,seq in enumerate(zip(*sequences), 1):
            if seq[0] == '-':
                i_length += 1
                i_offset += 1
            else:
                i_length = 0
            if i_length:
                num = f"{i-i_offset}_{i_length:03}"  
            else:
                num = f"{i-i_offset}"
            df[num] = [seq[0]] + ['=' if (s == seq[0]) else s for s in seq[1:]]
        return df

    df = compare_AAs(AA_alignment)
        #   Generate the mutations strings
    out = []
    for n,s in df.tail(-1).iterrows():
        out.append(f" {FORMATTING['G_FG']}{n                      :<{SEQ_L-1}}{FORMATTING['OFF']}|" +
                   f" {FORMATTING['R_FG']}{classify(classifier, s):<{VAR_L-1}}{FORMATTING['OFF']}|" +
                   f" {FORMATTING['Y_FG']}{mutations_string(df.iloc[0], s)   }{FORMATTING['OFF']}")
    return '\n'.join(out), df


def delta_convention_enforcer(lst: List[str]) -> List[str]:
    """Resolve a visualisation issue with the delta variant AA alignment"""
    if "E156G" in lst:
        i = lst.index("E156G")
        if lst[i + 1] == "∆157-158" and not match(".159", lst[i + 2]):
            lst[i:i + 2] = ["∆156-157", "R158G"]
    return lst


def load_reference(ref_file:List[str]) -> Tuple[SeqRecord, Seq]:
    """Load the reference DNA and AA sequences"""
    if not ref_file:
        raise FileNotFoundError("ERROR! Can't load reference sequence.")
    ref_DNA = read_fasta(ref_file[0])[0]
    ref_AA = translate_seq(ref_DNA.seq)
    return ref_DNA, ref_AA


def load_classifier(cls_file:List[str], train_file:List[str], ref:Seq) -> Optional[Classifier]:
    """Load and train the classifier"""
    if not cls_file:
        cls = None
    else:
        silent_print(f"Loading classifier from file '{cls_file[0]}'.")
        with open(cls_file[0], "rb") as f:
            cls = load(f)
                #   If an hardcoded threshold values is given, override the imported one
            cls.threshold = THRESHOLD or cls.threshold
    if train_file:
            #   Check whether the classifier is given, else instantiate a new one before training, then save it
        if not cls:
            silent_print(f"Initialising new classifier instance.")
            cls = Classifier(cls=RandomForestClassifier(**{'bootstrap'        : True
                                                          ,'max_depth'        : 80
                                                          ,'max_features'     : 'sqrt'
                                                          ,'min_samples_leaf' : 2
                                                          ,'min_samples_split': 5
                                                          ,'n_estimators'     : 20
                                                          ,'class_weight'     : 'balanced_subsample'}), ref=ref, threshold=THRESHOLD or 0.8)
        silent_print(f"Training classifier with data found in file '{train_file[0]}'.")
        cls.train(cls.readCSV(train_file[0]))
        cls.save((cls_out := f"{''.join(train_file[0].split('.')[:-1])}_classifier.pkl"))
        silent_print(f"Saving trained classifier as file '{cls_out}'.")
    return cls


def pipeline(*seq:SeqRecord, reference:SeqRecord, classifier:Optional[Classifier]) -> Tuple[str, DataFrame]:
    """The s_portein_typer core workflow, set up for parallelisation"""
    alignment1 = align_sequences(reference, *seq, reference=reference)
    translate_DNA_alignment(alignment1)
    alignment2 = align_sequences(*alignment1, reference=reference)
    return get_differences(alignment2, classifier)
    

def write_dataframe(file_name:str, df:DataFrame) -> None:
    """Write the given mutations DataFrame to a .csv file"""
    if file_name:
        silent_print(f"Saving alignment CSV file as '{file_name}'")
        strs_to_float = lambda ss: [float(s.replace('_', '.')) for s in ss]
        unique_indices = ~(df.duplicated() & df.index.duplicated())
        df[unique_indices].T.sort_index(key=strs_to_float).to_csv(file_name)


def main() -> None:
    """s_protein_typer main function workflow"""
    parser = ArgumentParser(description="Compares Sars-CoV2 spike protein sequences")
    parser.add_argument("-a", "--alignment",        nargs='?',            default="",                         metavar="alignment.csv",\
                            help="read a pre-existent alignment file, either .fasta or .csv, calculate and print its classification")
    parser.add_argument("-c", "--classifier",       nargs='?',            default="reference/classifier.pkl", metavar="classifier.pkl",\
                            help="relative or absolute path to the classifier binary .pkl file")
    parser.add_argument("-m", "--machine_readable", action="store_true",  default=False,
                            help="disable terminal formatting to produce a machine readable output")
    parser.add_argument("-o", "--output",           nargs='?',            default="",                         metavar="alignment.csv",\
                            help="after calculating the differences between the input sequences, save them to a separate .csv file")
    parser.add_argument("-s", "--slow",             action="store_true",  default=False,
                            help="perform a more accurate but slower Multiple Sequences Alignment")
    parser.add_argument("-r", "--reference",        nargs='?',            default="reference/ref_DNA.fasta",  metavar="ref_DNA.fasta",\
                            help="relative or absolute path to the reference .fasta file")
    parser.add_argument("-t", "--retrain",          nargs='?',            default="",                         metavar="variants.csv",\
                            help="relative or absolute path to the variants .csv database, to retrain the classifier")
    parser.add_argument("-v", "--verbose",          action="store_true",  default=False,
                            help="verbose mode")
    parser.add_argument("sequences",                nargs='*',            default=["data/*.fasta"],
                            help="relative or absolute path to the DNA sequences .fasta files")
    args = parser.parse_args()
        
    global VERBOSE, FORMATTING
    VERBOSE = args.verbose
    if args.machine_readable:
        silent_print(f"Deactivating terminal output formatting.")
        FORMATTING = defaultdict(str)
    
    ref_DNA, ref_AA = load_reference(glob(args.reference))
    classifier      = load_classifier(glob(args.classifier), glob(args.retrain), ref_AA)

    if not (alignment_file := glob(args.alignment)):
        if not (files := [g for f in args.sequences for g in glob(f)]):
            raise FileNotFoundError(f"ERROR! No file found for '{' '.join(args.sequences)}'.")
        sequences = read_fasta(*files, optimise=True)
            #   If the --slow flag is given, disable the parallelisation pipeline and run a single MSA
        if args.slow:
            silent_print(f"Performing multiple sequence alignment.")
            diffs, df = pipeline(*sequences, reference=ref_DNA, classifier=classifier)
        else:
            silent_print(f"Performing parallel sequence alignments.")
            diffss, dfs = zip(*Parallel(n_jobs=-1)(delayed(pipeline)(s, reference=ref_DNA, classifier=classifier) for s in sequences))
            diffs,  df  = '\n'.join(diffss), concat(dfs).fillna('-')
        write_dataframe(args.output, df)
    else:
        input_alignment = alignment_file[0]
        #   If an alignment file is given, try and read it as a .fasta or a .csv one
        try:
            alignment2 = list(read(input_alignment, "fasta"))
        except ValueError:
            alignment2 = read_csv(input_alignment, index_col=0).T
            alignment2.columns = alignment2.columns.astype("str")
        silent_print(f"Alignment found in '{input_alignment}'; computing differences.")        
        diffs, _ = get_differences(alignment2, classifier)

    silent_print(f"\r   \n{FORMATTING['SEQ_C']}{'Sequence ID':^{SEQ_L}}{FORMATTING['OFF']}"\
                +      f"|{FORMATTING['VAR_C']}{'Variant'    :^{VAR_L}}{FORMATTING['OFF']}"\
                +      f"|{FORMATTING['MUT_C']}{'Mutations'  :^15}{     FORMATTING['OFF']}")
    print(diffs)    
    silent_print(f"""with:""")
    silent_print(f"""    \r{FORMATTING['G_FG']}{'∆sequence' :>16}{FORMATTING['OFF']}: sequence derived from trimming
                         \r{FORMATTING['G_FG']}{'-sequence' :>16}{FORMATTING['OFF']}: sequence is too short to yield any meaningful result""")
    if classifier and not args.machine_readable:
        silent_print(f"""\r{FORMATTING['R_FG']}{'NA'        :>16}{FORMATTING['OFF']}: no class above the 0.5 confidence threshold
       \r{FORMATTING['R_FG']+FORMATTING['ITA']}{'class'     :>16}{FORMATTING['OFF']}: class confidence between the 0.5 and {classifier.threshold} thresholds
       \r{FORMATTING['R_FG']+FORMATTING['BOL']}{'class'     :>16}{FORMATTING['OFF']}: class confidence above the {classifier.threshold} threshold""")
    silent_print(f"""    \r{FORMATTING['Y_FG']}{'B123U'     :>15}{FORMATTING['OFF']} - mutation of aminoacid B into U at position 123
                         \r{FORMATTING['Y_FG']}{'∆123-456'  :>15}{FORMATTING['OFF']} - deletions from position 123 to position 456
                         \r{FORMATTING['Y_FG']}{'*123-456'  :>15}{FORMATTING['OFF']} - unrecognised residues from position 123 to position 456
                         \r{FORMATTING['Y_FG']}{'+123_BJU'  :>15}{FORMATTING['OFF']} - insertion of sequence BJU at position 123""")

if __name__ == "__main__":
    main()