def mfe(list sequences, str constraints, double temp , str dangles, int basepair=0):
    '''
    Calculate the MFE of a sequence using ViennaRNA as a module
    '''

    cdef object params = get_paramater_object(vienna_constants.material, temp, dangles, noLP = 1)
    cdef str seq_string = "&".join(sequences)

    cdef object rna = RNA.fold_compound(seq_string, params)

    if constraints:
        rna.constraints_add(constraints)

    cdef list mfe = rna.mfe()
    cdef str bracket_string = mfe[0]
    cdef double energy = round(mfe[1], 2)

    bp_x, bp_y = convert_bracket_to_numbered_pairs(bracket_string)

    cdef array.array mfe_basepairing_x = array.array('i', bp_x)
    cdef array.array mfe_basepairing_y = array.array('i', bp_y)
    cdef double mfe_energy = energy

    if basepair == 32012:
        print(sequences)

    return mfe_basepairing_x, mfe_basepairing_y, mfe_energy


