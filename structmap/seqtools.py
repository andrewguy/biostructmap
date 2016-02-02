from Bio import AlignIO

def _sliding_window(seq_align,window,isfile=True, fasta_out=False):
    '''
    Generator function that takes a multiple sequence alignment, and generates a
    Multiple Sequence Alignment over a sliding window.
    Input is either filehandle, or Bio.AlignIO multiple sequence alignment object.
    Set isfile=False if dealing with Bio.AlignIO MSA object.
    Output is either AlignIO object or fasta string.
    If fasta_out==False, then output will be AlignIO object.
    '''
    if isfile:
        alignments = AlignIO.read(seq_align,'fasta')
    else:
        alignments = seq_align
    #Length of alignments
    length = len(alignments[0])
    for i in range(length-window):
        alignment = alignments[:,i:i+window]
        if fasta_out:
            alignment = alignment.format('fasta')
        yield alignment

def _sliding_window_var_sites(seq_align,window,step=3,isfile=True,exclude_missing=True):
    '''
    Generator function that takes a multiple sequence alignment,
    and generates a Multiple Sequence Alignment over a sliding window, only
    including polymorphic sites in the alignment.
    Output is Bio.AlignIO alignment object.
    '''
    if isfile:
        alignments = AlignIO.read(seq_align,'fasta')
    else:
        alignments = seq_align
    #Length of alignments
    length = len(alignments[0])

    align_dict = _var_site(alignments)

    #Create first window
    initial_sites = {key:value for (key,value) in align_dict.items() if key<window}
    #Small hack to set type of 'initial_sites' variable if no alignments fall
    #within initial window
    initial_sites[-1] = alignments[:,0:0]

    alignment = _join_alignments(initial_sites)
    yield alignment
    #Add/remove sites from the end/start of window as appropriate.
    for i in range((length-window)//step):
        for j in range(step):
            if i + j in align_dict:
                alignment = alignment[:,1:]
            if i + j + window in align_dict:
                alignment = alignment + align_dict[i+j+window]
        yield alignment


def _var_site(alignment):
    '''
    Take a multiple sequence alignment object and return polymorphic sites in a
    dictionary object.
    Use this function to simplify the input to a tajima's D calculation.
    '''
    result = {}
    for i in range(len(alignment[0])):
        x = alignment[:,i]
        #Check if string contains a single character. Most efficient method found so far.
        if x != len(x) * x[0]:
            result[i] = alignment[:,i:i+1]
    return result

def _join_alignments(align_dict):
    '''
    Take a dictionary of multiple sequence alignments, and join according to
    dictionary key order (generally position in sequence).
    '''
    output = None
    for key in sorted(align_dict):
        if not output:
            output = align_dict[key]
        else:
            output = output + align_dict[key]
    return output
