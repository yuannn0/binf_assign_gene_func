from Bio.Align import substitution_matrices

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    #create table
    scores = []
    for i in range(seq1_len + 1):
        row = [0] * (seq2_len + 1)
        scores.append(row)

    #score first col
    for i in range(1, seq1_len + 1):
        scores[i][0] = scores[i - 1][0] + scoring_function(seq1[i - 1], "-")
    #score first row    
    for j in range(1, seq2_len + 1):
        scores[0][j] = scores[0][j - 1] + scoring_function("-", seq2[j - 1])

    #score the whole table
    for i in range(1, seq1_len + 1):
        for j in range(1, seq2_len + 1):
            diagonal = scores[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = scores[i - 1][j] + scoring_function(seq1[i - 1], "-")
            right = scores[i][j - 1] + scoring_function("-", seq2[j - 1])
            scores[i][j] = max(diagonal, up, right)
        
    #Traceback
    align_seq1 = []
    align_seq2 = []
    row = seq1_len
    col = seq2_len
    while row > 0 or col > 0:
        if row > 0 and col > 0 and scores[row][col] == scores[row-1][col-1] + scoring_function(seq1[row-1],seq2[col-1]):
            align_seq1.append(seq1[row - 1])
            align_seq2.append(seq2[col - 1])
            row -= 1
            col -= 1
        elif row > 0 and scores[row][col] == scores[row - 1][col] + scoring_function(seq1[row - 1], "-"):
            align_seq1.append(seq1[row - 1])
            align_seq2.append("-")
            row -= 1
        elif col > 0 and scores[row][col] == scores[row][col - 1] + scoring_function("-", seq2[col - 1]):
            align_seq2.append(seq2[col - 1])
            align_seq1.append("-")
            col -= 1

    align_seq1 = "".join(align_seq1[::-1])
    align_seq2 = "".join(align_seq2[::-1])

    return align_seq1, align_seq2, float(scores[seq1_len][seq2_len])

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    # create table
    scores = []
    for i in range(seq1_len + 1):
        table = [0] * (seq2_len + 1)
        scores.append(table)
        
    max_scores = 0
    position = (0,0)
    #score the whole table
    for i in range(1, seq1_len + 1):
        for j in range(1, seq2_len + 1):
            diagonal = scores[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = scores[i - 1][j] + scoring_function(seq1[i - 1], "-")
            right = scores[i][j - 1] + scoring_function("-", seq2[j - 1])
            scores[i][j] = max(0, diagonal, up, right)
            #track the highest score
            if scores[i][j] > max_scores:
                max_scores = scores[i][j]
                position = (i,j)

    #Traceback
    align_seq1 = []
    align_seq2 = []
    row, col = position
    while scores[row][col] > 0:
        if row > 0 and col > 0 and scores[row][col] == scores[row-1][col-1] + scoring_function(seq1[row-1],seq2[col-1]):
            align_seq1.append(seq1[row - 1])
            align_seq2.append(seq2[col - 1])
            row -= 1
            col -= 1
        elif row > 0 and scores[row][col] == scores[row - 1][col] + scoring_function(seq1[row - 1], "-"):
            align_seq1.append(seq1[row - 1])
            align_seq2.append("-")
            row -= 1
        elif col > 0 and scores[row][col] == scores[row][col - 1] + scoring_function("-", seq2[col - 1]):
            align_seq2.append(seq2[col - 1])
            align_seq1.append("-")
            col -= 1

    align_seq1 = "".join(align_seq1[::-1])
    align_seq2 = "".join(align_seq2[::-1])

    return align_seq1, align_seq2, float(max_scores)

## blosum62 scoring matrix 
blosum62 = substitution_matrices.load("BLOSUM62")

def scoring_function(aa_i,aa_j):
    #gap penalty
    if aa_i == "-" or aa_j == "-":
        return -5

    return blosum62[aa_i, aa_j]
