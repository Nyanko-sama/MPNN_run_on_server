import os 
import subprocess
import argparse 
from pathlib import Path
from typing import Optional
import tempfile
from Bio.PDB import MMCIFParser, PDBIO


parser = argparse.ArgumentParser()
parser.add_argument("--test_run", type=bool, default=False)
parser.add_argument("--input_folder", type=str, default="protein_conformations")
parser.add_argument("--output_folder", type=str, default="out")
parser.add_argument("--temps", type=list, default=[0.1, 0.2, 0.3])
parser.add_argument("--n_designs", type=int, default=1000)

args = parser.parse_args() 

def convert_cif_to_pdb(cif_path: Path, pdb_path: Optional[Path] = None, chain_id: Optional[str] = None) -> Path:
    """
    Convert CIF file to PDB format, optionally extracting only a specific chain.
    
    Args:
        cif_path: Path to the input CIF file
        pdb_path: Optional output PDB path. If None, creates a temp file.
        chain_id: Optional chain ID to extract (if None, extracts all chains)
    
    Returns:
        Path to the converted PDB file
    """
    if pdb_path is None:
        # Create a temporary PDB file with the same name as the cif file
        pdb_path = Path(tempfile.mkstemp(suffix='.pdb', prefix=f'{cif_path.stem}_')[1])
    
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('structure', str(cif_path))
        
        io = PDBIO()
        if chain_id:
            # Extract only the specified chain
            class ChainSelector:
                def __init__(self, chain_id):
                    self.chain_id = chain_id
                def accept_model(self, model):
                    return 1
                def accept_chain(self, chain):
                    return chain.id == self.chain_id
                def accept_residue(self, residue):
                    return 1
                def accept_atom(self, atom):
                    return 1
            
            selector = ChainSelector(chain_id)
            io.set_structure(structure)
            io.save(str(pdb_path), select=selector)
        else:
            io.set_structure(structure)
            io.save(str(pdb_path))
        
        return pdb_path
    except Exception as e:
        raise RuntimeError(f"Failed to convert CIF to PDB: {e}") from e

def get_pdb_path_for_mpnn(input_path: Path, temp_dir: Optional[Path] = None) -> Path:
    """
    Get a PDB file path for ProteinMPNN. Converts CIF to PDB if needed.
    
    Args:
        input_path: Path to input file (PDB or CIF)
        temp_dir: Optional directory for temporary PDB files
    
    Returns:
        Path to PDB file (original if input is PDB, converted if input is CIF)
    """
    if input_path.suffix.lower() == '.pdb':
        return input_path
    elif input_path.suffix.lower() == '.cif':
        if temp_dir:
            temp_dir.mkdir(parents=True, exist_ok=True)
            pdb_path = temp_dir / f"{input_path.stem}.pdb"
        else:
            pdb_path = None  # Will create temp file
        #get chain id from cif file name, example: cndt_1ib1E_aligned_to_1kuvA_0.fa -> protein and chain id is 1ib1E -> chain id is E 
        chain_id = input_path.stem.split("_")[1][-1] if len(input_path.stem.split("_")) > 1 else None
        print(f"Chain id: {chain_id}")
        return convert_cif_to_pdb(input_path, pdb_path, chain_id=chain_id)
    else:
        raise ValueError(f"Unsupported file format: {input_path.suffix}. Expected .pdb or .cif")

def run_proteinmpnn_with_auto_convert(
    input_path: Path,
    output_dir: Path,
    num_seq: int = 10,
    temp: str = "0.1",
    temp_dir: Optional[Path] = None,
    cleanup: bool = True
) -> bool:
    """
    Run ProteinMPNN with automatic CIF to PDB conversion.
    
    Args:
        input_path: Path to input file (PDB or CIF)
        output_dir: Output directory for ProteinMPNN
        num_seq: Number of sequences to generate
        temp: Sampling temperature
        temp_dir: Directory for temporary PDB files (if None, uses system temp)
        cleanup: Whether to clean up temporary PDB files after running
    
    Returns:
        True if successful, False otherwise
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert CIF to PDB if needed
    pdb_path = get_pdb_path_for_mpnn(input_path, temp_dir)
    is_temp_file = (pdb_path != input_path)
    
    try:
        # Build ProteinMPNN command
        cmd = [
            'python', 'ProteinMPNN/protein_mpnn_run.py',
            '--pdb_path', str(pdb_path),
            '--out_folder', str(output_dir),
            '--num_seq_per_target', str(num_seq),
            '--sampling_temp', temp,
        ]
        
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR running ProteinMPNN: {e.stderr}")
        return False
    finally:
        # Clean up temporary PDB file if we created one
        if cleanup and is_temp_file and pdb_path.exists():
            try:
                pdb_path.unlink()
            except Exception:
                pass



# Main 
def main():
    args = parser.parse_args()
    if args.test_run:
        print("Running in test mode")
        # run only for 1 conformation from protein_conformations folder (the conformation is in subfolder)
        for protein_folder in os.listdir(args.input_folder):
            if os.path.isdir(os.path.join(args.input_folder, protein_folder)):
                for conformation_folder in os.listdir(os.path.join(args.input_folder, protein_folder)):
                    if os.path.isdir(os.path.join(args.input_folder, protein_folder, conformation_folder)):
                        for temp in args.temps: 
                            print(f"Running for {protein_folder} {conformation_folder} {temp}")
                            input_path = Path(os.path.join(args.input_folder, protein_folder, conformation_folder))
                            output_path = Path(os.path.join(args.output_folder, protein_folder, conformation_folder, f"temp_{temp}"))
                            run_proteinmpnn_with_auto_convert(input_path, output_path, args.n_designs, str(temp))
                        break
                break

    else:
        print("Running in full mode")
        # run for all conformations in protein_conformations folder
        for protein_folder in os.listdir(args.input_folder):
            if os.path.isdir(os.path.join(args.input_folder, protein_folder)):
                for conformation_folder in os.listdir(os.path.join(args.input_folder, protein_folder)):
                    if os.path.isdir(os.path.join(args.input_folder, protein_folder, conformation_folder)):
                        for temp in args.temps: 
                            print(f"Running for {protein_folder} {conformation_folder} {temp}")
                            input_path = Path(os.path.join(args.input_folder, protein_folder, conformation_folder))
                            output_path = Path(os.path.join(args.output_folder, protein_folder, conformation_folder, f"temp_{temp}"))
                            
                            run_proteinmpnn_with_auto_convert(input_path, output_path, args.n_designs, str(temp))


if __name__ == "__main__":
    main()