from dataclasses import dataclass

@dataclass
class OSTIRResult():
    """Class for storing the results of an OSTIR calculation"""
    # These are the results that most people care about
    name: str = ""
    start_codon: str = ""
    start_position: int = 0
    expression: float = 0.0
    RBS_distance_bp: int = 0
    dG_total: float = 0.0
    dG_rRNA_mRNA: float = 0.0
    dG_mRNA: float = 0.0
    dG_spacing: float = 0.0
    dG_standby: float = 0.0
    dG_start_codon: float = 0.0

    # These are the results that are not as important as the above
    mRNA_structure: str = ""
    mRNA_rRNA_uncorrected_structure: str = ""
    mRNA_rRNA_corrected_structure: str = ""
    most_5p_mRNA: str = ""
    longest_helix: str = ""
    dG_start_energy: float = 0.0
    helical_loop: str = ""
    bulge_loop: str = ""
    min_bp_prob: float = 0.0
    kinetic_score: float = 0.0

    def __getitem__(self, key):
        if key in self.__annotations__:
            return getattr(self, key)
        raise KeyError(key)

    def results(self):
        """Returns a dictionary of the core results"""
        return {
            'name': self.name,
            'start_codon': self.start_codon,
            'start_position': self.start_position,
            'expression': self.expression,
            'RBS_distance_bp': self.RBS_distance_bp,
            'dG_total': self.dG_total,
            'dG_rRNA:mRNA': self.dG_rRNA_mRNA,
            'dG_mRNA': self.dG_mRNA,
            'dG_spacing': self.dG_spacing,
            'dG_standby': self.dG_standby,
            'dG_start_codon': self.dG_start_codon,
        }

    def addional_results(self):
        """Returns a dictionary additional results"""
        return {
            'mRNA_structure': self.mRNA_structure,
            'mRNA_rRNA_uncorrected_structure': self.mRNA_rRNA_uncorrected_structure,
            'mRNA_rRNA_corrected_structure': self.mRNA_rRNA_corrected_structure,
            'most_5p_mRNA': self.most_5p_mRNA,
            'longest_helix': self.longest_helix,
            'dG_start_energy': self.dG_start_energy,
            'helical_loop': self.helical_loop,
            'bulge_loop': self.bulge_loop,
            'min_bp_prob': self.min_bp_prob,
            'kinetic_score': self.kinetic_score,
        }

@dataclass
class ostir_task:
    sequence: str
    start: int
    end: int
    name: None
    rrna: str='ACCTCCTTA'
    decimals: int=2
    circular: bool=False
    constraints: str=None

