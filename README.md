# Euler-Simulation
# Topological Analysis with Random Binary Matrices

This project provides a Python implementation for generating and analyzing the topological properties of binary matrices using simplicial complexes and Betti numbers. The code includes a simulation to calculate Euler characteristics and ensures consistency between Betti numbers and simplicial counts.

## Features

- Generates random binary matrices of specified dimensions.
- Identifies simplices (points, edges, triangles, tetrahedra) based on the structure of black pixels in the binary matrix.
- Constructs sparse boundary matrices for topological analysis.
- Calculates Betti numbers (β₀, β₁, β₂, β₃) and Euler characteristics.
- Validates Euler characteristic consistency through simulations.

## Requirements

To run this project, ensure the following Python libraries are installed:

- `numpy`
- `scipy`
- `sympy`

You can install these libraries using pip:

```bash
pip install numpy scipy sympy
Usage
Running the Simulation
Clone this repository or download the source code.
Run the Python script using a terminal or IDE.
Input the number of iterations for the simulation when prompted.
bash
Kodu kopyala
python topological_analysis.py
Example Output
The program generates random binary matrices, calculates topological invariants, and validates Euler consistency. It outputs the following:

Betti numbers for each dimension.
Euler characteristic values (from Betti numbers and simplices).
A consistency check for the Euler characteristic.
If all iterations are consistent, the program prints:

sql
Kodu kopyala
Euler Characteristic Consistency is TRUE for all iterations!
Total execution time: XX.XX seconds
Simulation complete.
Input Parameters
Matrix Size: Fixed at 30x30.
Black Pixel Count Range: Randomly chosen between 200 and 400 pixels.
Iterations: Defined by user input.
Code Structure
find_black_pixel_neighbors: Finds neighbors of black pixels in a binary matrix.
generate_simplices: Identifies 0-simplices, 1-simplices, 2-simplices, and 3-simplices.
generate_boundary_matrices_sparse: Constructs sparse boundary matrices for simplices.
approximate_smith_normal_form: Approximates the Smith Normal Form of a matrix.
calculate_betti_numbers_and_euler_characteristic: Computes Betti numbers and Euler characteristic, and validates consistency.
generate_random_binary_matrix: Creates random binary matrices with specified black pixel counts.
Example Simulation
The program simulates the process on a random binary matrix and validates the consistency of Euler characteristics across multiple iterations.

python
Kodu kopyala
Enter the number of iterations: 100
Euler Characteristic Consistency is TRUE for all iterations!
Total execution time: 45.32 seconds
Simulation complete.
Contribution
Contributions are welcome! If you'd like to improve or extend the functionality, please submit a pull request.

License
This project is licensed under the MIT License. See the LICENSE file for more information.
