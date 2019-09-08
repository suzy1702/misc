import Parsers
import copy as cp

class params(Parsers.parse_input):
    def simulation(self):
        self.steps = self.num_steps//self.stride
        self.steps_per_ensemble = self.steps//self.num_ensembles
    def lattice(self):
        lattice = Parsers.parse_lattice(self)
        self.num_atoms = lattice.num_atoms
        self.num_unitcells = lattice.num_unitcell
        self.num_basis = lattice.num_basis
        self.atom_ids = lattice.atom_ids
        self.unit_cells = lattice.unit_cells
        self.basis_pos = lattice.basis_pos
        self.masses = lattice.masses

