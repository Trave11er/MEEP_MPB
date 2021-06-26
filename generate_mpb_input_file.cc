// Generates .ctl (and .cell file to quickly visualize with Jmol) for a hexagon/parallelogram/triangle inside trivial matrix
// Think of the supercell as single cells (with cylinders or vacuum); single cell vectors are (sqrt(3)/2, 1/2) and (sqrt(3)/2, -1/2)
#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cstdlib>
using namespace std;

// run all the tests; have also checked that geometries for different numbers of top., trivial and vacuum cells are correct and reproduced result for bulk crystal
void unitary_test(double topological_cluster_radius, double trivial_cluster_radius, int supercell_size, double topological_radius) {
    assert(0.0 < trivial_cluster_radius && trivial_cluster_radius < 1.0/3); 
    assert(1.0/3 < topological_cluster_radius && topological_cluster_radius < 1.0/2); 
    assert(supercell_size > 2*topological_radius+0.1);
    assert((supercell_size-1)%2==0);
}

// occupy a single unit cell of material
void occupy_cell(double cluster_radius, int material_cell_num, int lat_vec_one_num, int lat_vec_two_num, double one_cylinder_positions[], double two_cylinder_positions[]){
    one_cylinder_positions[6*material_cell_num+0]=double(lat_vec_one_num)+(1.0/3.0)+(1.0/3.0);
    one_cylinder_positions[6*material_cell_num+1]=double(lat_vec_one_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    one_cylinder_positions[6*material_cell_num+2]=double(lat_vec_one_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    one_cylinder_positions[6*material_cell_num+3]=double(lat_vec_one_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    one_cylinder_positions[6*material_cell_num+4]=double(lat_vec_one_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    one_cylinder_positions[6*material_cell_num+5]=double(lat_vec_one_num)+(1.0/3.0)+(1.0/3.0);

    two_cylinder_positions[6*material_cell_num+0]=double(lat_vec_two_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    two_cylinder_positions[6*material_cell_num+1]=double(lat_vec_two_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    two_cylinder_positions[6*material_cell_num+2]=double(lat_vec_two_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    two_cylinder_positions[6*material_cell_num+3]=double(lat_vec_two_num)+(1.0/3.0)+(1.0/3.0);
    two_cylinder_positions[6*material_cell_num+4]=double(lat_vec_two_num)+(1.0/3.0)+(1.0/3.0);
    two_cylinder_positions[6*material_cell_num+5]=double(lat_vec_two_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
}

int main() {
    // Parameters for geometry
    // From Wu, Hu et al 1/2.9, 1/3.125 and epsilon(si)=11.7
    double topological_cluster_radius = 1.0/2.9; // radius of topological 6-cylinder cluster; not to be confused with cylinder radius
    double trivial_cluster_radius =1.0/3.125; // radius of trivial 6-cylinder cluster; not to be confused with cylinder radius
    int supercell_size=21; // size of supecell fully covereged with trivial clusters; should be odd
    int topological_radius=5; // # of top. cells in rad. dir. (along a side) of embedded topological supercluster at the origin; 6 cylinder each
    char particle_type='h'; // # h for hexagon (p for parallelogram and t for triangle) and r for rectangle (has hard-coded size)

    // Parameters for MPB calculation
    int number_of_bands = 1350;

    unitary_test(topological_cluster_radius,trivial_cluster_radius, supercell_size, topological_radius);
    double delta = 0.01; // small # for numerical comparissons

    // Initialise data structures
    int material_cells = supercell_size*supercell_size;
    int dummy = 0;
    int topological_count=0;
    double vector_norm_sq = 0;
    double x_eff, y_eff, phi, factor, eff_radius; // variables for chiral calculation

    // Initialise more arrays for cylinder coordinates
    double * one_cylinder_positions;
    double * two_cylinder_positions;
    one_cylinder_positions = new double [6*material_cells];
    two_cylinder_positions = new double [6*material_cells];

    // Save .cell file - can upload into Jmol to check the geometry
    ofstream write_cell;
    write_cell.open("output_silicon_embedded.cell");
    write_cell << "%block lattice_cart \n";
    write_cell << 10*supercell_size*sqrt(3)/2 << " " << 10*supercell_size*0.5 << " 0 \n";
    write_cell << 10*supercell_size*sqrt(3)/2 << " " << -10*supercell_size*0.5 << " 0 \n";
    write_cell << "0 0 10 \n";
    write_cell << "%endblock lattice_cart \n %block positions_frac \n";

    // Output cylinder positions in .ctl file for use in MPB
    ofstream write_ctl;
    write_ctl.open("output_silicon_embedded.ctl");
    // Define epsilon here 11.7 in Wu et al
    write_ctl << "(define cylinder_epsilon 11.7) \n (set! num-bands " << number_of_bands << ") \n";
    write_ctl << "(set! geometry-lattice (make lattice (size " << supercell_size << " " << supercell_size << " no-size) \n (basis1 (/ (sqrt 3) 2) 0.5) (basis2 (/ (sqrt 3) 2) -0.5))) \n";
    write_ctl << "(define cylinder_radius " << topological_cluster_radius/3 << " ) \n";
    write_ctl << "(set! geometry (list \n";

    // Occupy the position arrays with values
    for (int i=0; i<supercell_size; i++) {
        for (int j=0; j<supercell_size; j++) {
            // topological superculsters is in the centre of the supercell; at (size*.5+.5,size*.5+.5)
            vector_norm_sq=(i-supercell_size*0.5+0.5)*(i-supercell_size*0.5+0.5)+(j-supercell_size*0.5+0.5)*(j-supercell_size*0.5+0.5)+(i-supercell_size*0.5+0.5)*(j-supercell_size*0.5+0.5);
            switch (particle_type){
                case 'h':
                    // If the cell is within the topological radius; 
                    if (vector_norm_sq > topological_radius*topological_radius+delta ) {
                        occupy_cell(trivial_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        dummy++;
                    } else  {
                       if ((j==13 && i==12) || (j==12 && i==13) || (j==12 && i==12)) { 
                         occupy_cell(topological_cluster_radius+1e-5, dummy, i, j, one_cylinder_positions, two_cylinder_positions); 
                        } else {
                            occupy_cell(topological_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        }
                    topological_count++;
                    dummy++;
                    }   // If the cell is between the topological and material radius, i.e. it is trivial
                    break;
                        case 't':
                    if (i+j <= topological_radius) {
                        occupy_cell(topological_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        topological_count++;
                        dummy++;
                    }   // If the cell is between the topological and material radius, i.e. it is trivial
                    else {
                        occupy_cell(trivial_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        dummy++;
                    }
                    break;
                        case 'p':
                    if (i <= topological_radius && j <= topological_radius) {
                        occupy_cell(topological_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        topological_count++;
                        dummy++;
                    }   // If the cell is between the topological and material radius, i.e. it is trivial
                    else {
                        occupy_cell(trivial_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        dummy++;
                    }
                    break;
                        case 'r':
                    if ( ((sqrt(3)*0.5*(i+j))<24) && ((sqrt(3)*0.5*(i+j))>4) && ((0.5*(i-j))<1.5) && ((0.5*(i-j))>-1) ) {
                        occupy_cell(topological_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        topological_count++;
                        dummy++;
                    }   // If the cell is between the topological and material radius, i.e. it is trivial
                    else {
                        occupy_cell(trivial_cluster_radius, dummy, i, j, one_cylinder_positions, two_cylinder_positions);
                        dummy++;
                    }
                    break;
                        default:
                    cout << "Wrong particle type \n";
                    }
                    for (int k=0; k<6; k++) {
                        // for .cell file 
                        write_cell << "C " << one_cylinder_positions[6*dummy-6+k]/supercell_size << " " << two_cylinder_positions[6*dummy-6+k]/supercell_size << " 0 \n";
                        // for .ctl file
                        write_ctl << "(make cylinder \n (center " << one_cylinder_positions[6*dummy-6+k]-31.0/3 << " " << two_cylinder_positions[6*dummy-6+k]-31.0/3 << " 0) (radius cylinder_radius) (height infinity) \n (material (make dielectric (epsilon cylinder_epsilon)))) \n";
                    }
            }
        }

        assert(dummy==material_cells);
        cout << topological_radius << " is the radius of topological supercluster of shape " << particle_type << " with constituent clusters of radius " << topological_cluster_radius << endl <<
            "trivial cells have radius " << trivial_cluster_radius << endl
            << supercell_size << " is the supercell size and total number of clusters is " << material_cells << " of which topological are " << 
            topological_count << endl;
        write_cell << "%endblock positions_frac \n";
        write_cell.close();
        write_ctl << ")) \n";
        write_ctl << "(set! k-points (list (vector3 0 0 0) )) \n (set! resolution 32) \n (run-tm (output-at-kpoint (vector3 0 0 0) fix-efield-phase output-efield-z))";
        write_ctl.close();

        delete[] one_cylinder_positions;
        delete[] two_cylinder_positions;
        return 0;
    }

