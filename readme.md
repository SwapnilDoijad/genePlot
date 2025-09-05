# GenePlot

`GenePlot` is a lightweight bioinformatics pipeline built on [gggenomes](https://thackl.github.io/gggenomes/) for **comparative visualization of genes** across multiple genomic segments. It automates the extraction of coding sequences from GFF files, aligns homologous genes, and generates **comparative gene plots** to highlight synteny and structural variations. 

The order of genomic regions is determined by **Shared Genomic Content (SGC)**, calculated as:

> **SGC = ANI × shared region / 2**

---

## Key Features
- **Gene Synteny Visualization**: Compare gene order and structure across genomes.
- **Automated Workflow**: Extract coding sequences, align genes, and generate plots with minimal manual intervention.
- **Customizable Outputs**: Generate plots with original gene order or reordered by non-duplicated gene names to highlight missing genes.

---

## Example Plots

### Gene Synteny Plot (Original Gene Order)
The shaded areas between genomes represent pairwise BLAST matches.

<div align="center">
  <img src="test_results/genePlot.coordinates.png" alt="Gene Synteny Plot (Original Order)" width="700" height="600">
</div>

---

### Gene Synteny Plot (Reordered by Non-Duplicated Names)
This plot highlights missing genes (white regions).

<div align="center">
  <img src="test_results/genePlot.genes.png" alt="Gene Synteny Plot (Reordered)" width="700" height="600">
</div>

---

## Installation

Follow these steps to install and set up `GenePlot`:

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-repo/genePlot.git
   ```

2. **Navigate to the Project Directory**:
   ```bash
   cd genePlot
   ```

3. **Make Scripts Executable**:
   ```bash
   chmod +x installation.sh
   chmod +x ./ContigClass.sh
   ```

4. **Set Up the Conda Environment**:
   ```bash
   conda env create -f scripts/environment.yml
   conda activate genePlot
   ```

---

## Usage

Run the `genePlot.sh` script with the following options:

```bash
conda activate genePlot
./genePlot.sh --help
```

### Options:
- `-d, --data_dir <directory>`: Input data directory containing GFF files.
- `-o, --output_dir <directory>`: Output directory for results.
- `-c, --cpus <number>`: Number of CPUs to use (optional, default: 4).

### Example:
```bash
./genePlot.sh -d /path/to/gff/files -o /path/to/output -c 4
```

---

## Results

After running the script, the following files will be generated in the specified output directory:

- **`genePlot.coordinates.png`**: Gene synteny plot with original gene order.
- **`genePlot.genes.png`**: Gene synteny plot reordered by non-duplicated names.
- Additional intermediate files and logs for debugging.

---

## Author

🧑‍💻 **Swapnil Doijad**  
📧 Email: [swapnil.doijad@gmail.com](mailto:swapnil.doijad@gmail.com)

---

## Support

If you encounter bugs or have feature requests, please open an issue in the repository. Contributions are welcome!