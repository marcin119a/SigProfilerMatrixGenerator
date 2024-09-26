import itertools

def perm(n, seq):
    """
    Generates a list of all available permutations of n-mers.

    Parameters:
                       n  -> length of the desired permutation string
                     seq  -> list of all possible string values

    Returns:
              permus  -> list of all available permutations
    """
    permus = []
    for p in itertools.product(seq, repeat=n):
        permus.append("".join(p))
    return permus

size = 5
bases = ["A", "C", "G", "T"]
mut_types_initial = perm(size, "ACGT")
mut_types = []
tsb = ["T", "U", "N", "B"]
for tsbs in tsb:
        for mut in mut_types_initial:
            current_base = mut[int(size / 2)]
            if current_base == "C" or current_base == "T":
                for base in bases:
                    if base != current_base:
                        mut_types.append(
                            tsbs
                            + ":"
                            + mut[0: int(size / 2)]
                            + "["
                            + current_base
                            + ">"
                            + base
                            + "]"
                            + mut[int(size / 2) + 1:]
                        )
mutation_types_tsb_context = []
mutation_types = [
    "CC>AA",
    "CC>AG",
    "CC>AT",
    "CC>GA",
    "CC>GG",
    "CC>GT",
    "CC>TA",
    "CC>TG",
    "CC>TT",
    "CT>AA",
    "CT>AC",
    "CT>AG",
    "CT>GA",
    "CT>GC",
    "CT>GG",
    "CT>TA",
    "CT>TC",
    "CT>TG",
    "TC>AA",
    "TC>AG",
    "TC>AT",
    "TC>CA",
    "TC>CG",
    "TC>CT",
    "TC>GA",
    "TC>GG",
    "TC>GT",
    "TT>AA",
    "TT>AC",
    "TT>AG",
    "TT>CA",
    "TT>CC",
    "TT>CG",
    "TT>GA",
    "TT>GC",
    "TT>GG",
]

for base in bases:
    for mut in mutation_types:
        for base2 in bases:
            for base3 in tsb:
                mutation_types_tsb_context.append(
                    "".join([base3, ":", base, "[", mut, "]", base2])
                )