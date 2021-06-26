// Generates .ctl and .cell file (to check with Jmol) for MEEP of either hexagonal topological
// particle in trivial matrix or two respective half planes
#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cstdlib>
using namespace std;

// run all the tests; have also checked that geometries for different numbers of top., trivial and vacuum cells are correct and reproduced result for bulk crystal
void unitary_test(double topological_cluster_radius, double trivial_cluster_radius, int supercell_size_x, int supercell_size_y,int  particle_radius) {
    assert(0.0 < trivial_cluster_radius && trivial_cluster_radius < 0.333334); 
    assert(0.333334 < topological_cluster_radius && topological_cluster_radius < 0.5); 
    assert (particle_radius < supercell_size_x/2 && particle_radius < supercell_size_y/2);
}

// occupy a single unit cell of material
void occupy_cell(double cluster_radius, int lat_vec_one_num, int lat_vec_two_num, double one_cylinder_positions[], double two_cylinder_positions[]){
    one_cylinder_positions[0]=double(lat_vec_one_num)+(1.0/3.0)+(1.0/3.0);
    one_cylinder_positions[1]=double(lat_vec_one_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    one_cylinder_positions[2]=double(lat_vec_one_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    one_cylinder_positions[3]=double(lat_vec_one_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    one_cylinder_positions[4]=double(lat_vec_one_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    one_cylinder_positions[5]=double(lat_vec_one_num)+(1.0/3.0)+(1.0/3.0);

    two_cylinder_positions[0]=double(lat_vec_two_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    two_cylinder_positions[1]=double(lat_vec_two_num)+(1.0/3.0-cluster_radius)+(1.0/3.0);
    two_cylinder_positions[2]=double(lat_vec_two_num)+(1.0/3.0+cluster_radius)+(1.0/3.0);
    two_cylinder_positions[3]=double(lat_vec_two_num)+(1.0/3.0)+(1.0/3.0);
    two_cylinder_positions[4]=double(lat_vec_two_num)+(1.0/3.0)+(1.0/3.0);
    two_cylinder_positions[5]=double(lat_vec_two_num)+(1.0/3.0+cluster_radius)+(1.0/3.0); 
}

int main() {
    // Parameters for geometry
    // From Wu, Hu et al 1/2.9, 1/3.125 and epsilon(Si)=11.7
    double topological_cluster_radius = 1.0/2.9; // radius of topological 6-cylinder cluster; not to be confused with cylinder radius
    double trivial_cluster_radius = 1.0/3.125; // radius of trivial 6-cylinder cluster; not to be confused with cylinder radius
    int dummy_size_x=60; // size of cell in x-direction
    int dummy_size_y=60; // size of cell in y-direction
    int supercell_size_x=13; // size of cell in x-direction
    int supercell_size_y=13; // size of cell in y-direction
    int particle_radius =3; // radius of the particle (only used for particle type 'h')
    char particle_type='h';  // 'h' for hexagonal, 'r' for half-plane
 
    double x,y;
    unitary_test(topological_cluster_radius,trivial_cluster_radius, supercell_size_x, supercell_size_y, particle_radius);
    double delta = 0.01; // small # for numerical comparissons
    double cell_file_shift = 0.0; // set to 0.5 to ease visualization with jmol (but has different coordinates from .ctl then)
    srand(0);

    // Initialise data structures
    int dummy = 0;
    int topological_count=0;

    // Initialise more arrays for cylinder coordinates
    double * one_cylinder_positions; // stores x-cpt of position
    double * two_cylinder_positions; // stores y-cpt of position
    one_cylinder_positions = new double [6];
    two_cylinder_positions = new double [6];

    // Save .cell file - can upload into Jmol to check the geometry
    ofstream write_cell;
    write_cell.open("output_meep_coords.cell");
    write_cell << "%block lattice_cart \n";
    write_cell << 10*supercell_size_x << " " << 0 << " 0 \n";
    write_cell << 0 << " " << 10*supercell_size_y << " 0 \n";
    write_cell << "0 0 10 \n";
    write_cell << "%endblock lattice_cart \n %block positions_frac \n";

    // Output cylinder positions in .ctl file for use in MPB
    ofstream write_ctl;
    write_ctl.open("output_meep_coords.ctl");
    write_ctl << "(define cylinder_epsilon 11.7) \n"; 
    write_ctl << "(set! geometry-lattice (make lattice (size " << supercell_size_x << " " << supercell_size_y << " no-size))) \n";
    write_ctl << "(define cylinder_radius " << topological_cluster_radius/3 << " ) \n";
    write_ctl << "(set! geometry (list \n";

    // Occupy the position arrays with values
    for (int i=-dummy_size_x/2+1; i<dummy_size_x/2-1; i++) {
        for (int j=-dummy_size_y/2+1; j<dummy_size_y-1; j++) {
            x=i+j*0.5;
            y=j*0.5*sqrt(3);
            switch(particle_type) {
            // To imitate topological ribbon
                case 'r':
            if ((x < supercell_size_x/2-1 && x > -supercell_size_x/2+1 && y < supercell_size_y/2-1 && y > -supercell_size_y/2+1)){
                // set y<0 to make lower half-plane topological 
                if (y < 0) {
                    occupy_cell(topological_cluster_radius, i, j, one_cylinder_positions, two_cylinder_positions);
                    topological_count++;
                    dummy++;
                } else  {
                    occupy_cell(trivial_cluster_radius, i, j, one_cylinder_positions, two_cylinder_positions);
                    dummy++;
                }
                for (int k=0; k<6; k++) {
                    // for .cell file
                    write_cell << "C " << (one_cylinder_positions[k]+0.5*two_cylinder_positions[k])/supercell_size_x << " " << 0.5*sqrt(3)*two_cylinder_positions[k]/supercell_size_y << " 0 \n";
                    // for .ctl file
                    write_ctl << "(make cylinder \n (center " << one_cylinder_positions[k]+0.5*two_cylinder_positions[k] << " " << 0.5*sqrt(3)*two_cylinder_positions[k] << " 0) (radius cylinder_radius) (height infinity) \n (material (make dielectric (epsilon cylinder_epsilon)))) \n";
                }
            }
            break;
            // For hexagonal topological particle
                case 'h':
            if (x < supercell_size_x/2-1 && x > -supercell_size_x/2+1 && y < supercell_size_y/2-1 && y > -supercell_size_y/2+1) {
                // upper half-plane is topological
                if (x*x+y*y > particle_radius*particle_radius) {
                    occupy_cell(trivial_cluster_radius, i, j, one_cylinder_positions, two_cylinder_positions);
                    dummy++;
                } 
                else  {
                    // add an air or trivial defect only works for particle radius = 5
                    occupy_cell(topological_cluster_radius, i, j, one_cylinder_positions, two_cylinder_positions);
                    topological_count++;
                    dummy++;
                }

                for (int k=0; k<6; k++) {
                    // for .cell file
                    write_cell << "C " << cell_file_shift+(one_cylinder_positions[k]+0.5*two_cylinder_positions[k])/supercell_size_x << " " << cell_file_shift+0.5*sqrt(3)*two_cylinder_positions[k]/supercell_size_y << " 0 \n";
                    // for .ctl file
                    write_ctl << "(make cylinder \n (center " << one_cylinder_positions[k]+0.5*two_cylinder_positions[k] << " " << 0.5*sqrt(3)*two_cylinder_positions[k] << " 0) (radius cylinder_radius) (height infinity) \n (material (make dielectric (epsilon cylinder_epsilon)))) \n";
                }
            }
            break;
            default:
            cout << "wront particle type" << endl;
            }
        }
    }
    cout <<  "topological supercluster of shape " << particle_type << " with constituent clusters of radius " << topological_cluster_radius << endl <<
        "trivial cells have radius " << trivial_cluster_radius << endl
        << supercell_size_x << "x" << supercell_size_y << " is the supercell size and number of topologicla clusters is "  << 
        topological_count << endl;
    write_cell << "%endblock positions_frac \n";
    write_cell.close();
    write_ctl << ")) \n"; // to close the geometry list
    write_ctl << "(set! pml-layers (list (make pml (thickness 1.0)))) \n " << 
        "(set! resolution 20) \n " << 
        "(define-param nfreq 200) \n " <<
        "(define-param fcen 0.473) \n " <<
        "(define-param df 0.05) \n " <<
        "(set! sources (list \n " << 
        "(make source \n " << 
        "(src (make gaussian-src (frequency fcen) (fwidth df))) \n "
        "(component Hx) \n " 
        "(center 0.5 0.0) \n " << 
        "(amplitude 1)) " << 
        "(make source \n " << 
        "(src (make gaussian-src (frequency fcen) (fwidth df))) \n "
        "(component Hy) \n "
        "(center 0.5 0.0) \n " << 
        "(amplitude (exp (* 0+1i 1.570796327)) )) " <<  // To make a circular dipole
        ")) \n " <<
        "(define left \n " << // first monitor
        "(add-flux fcen df nfreq \n "<<
        "(make flux-region \n " << 
        "(center -3 0 1) (size 0 2 2)))) \n \n " <<
        "(define right \n " << // second monitor
        "(add-flux fcen df nfreq \n "<<
        "(make flux-region \n " << 
        "(center 4 0 1) (size 0 2 2)))) \n \n " <<
        "(run-sources+ 500 " << 
        "(at-beginning output-epsilon)) \n \n " << 
 //       "(run-sources+ \n (stop-when-fields-decayed 50 Ez \n (vector3 0.5 4.6) 1e-3)) \n \n " <<
        "(display-fluxes left right)";
    write_ctl.close();

    delete[] one_cylinder_positions;
    delete[] two_cylinder_positions;
    return 0;
}
