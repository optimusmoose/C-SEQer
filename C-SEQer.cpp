/*
	This work is a derivative of the original SweetSEQer algorithm licensed to Oliver Serang (2012):
		https://pubmed.ncbi.nlm.nih.gov/23443135/
		
	Copyright 2021 Christopher Burgoyne
	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at
		http://www.apache.org/licenses/LICENSE-2.0
	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<iomanip>
#include<sstream>

using namespace std;

#define print_precision 6
#include "main_functions.h"
#include "MGF.h"
#include "DeNovo.h"

int main(int argc, char *args[]){
	if (argc == 9) {

		double error_tolerance;
        double peak_intensity_tolerance;
        int min_peptide_length;
        int min_glycan_length;
        double min_mz_for_glycan;
		bool print_output;
		bool draw_figures;
		string mgf_fname;
		ifstream mgf_file;
		
		capture_parameters(args, &error_tolerance, &peak_intensity_tolerance, &min_peptide_length, &min_glycan_length, &min_mz_for_glycan, &print_output, &draw_figures, &mgf_fname, mgf_file);
		print_parameters(error_tolerance, peak_intensity_tolerance, min_peptide_length, min_glycan_length, min_mz_for_glycan, draw_figures);
		
		// Counters
		int number_spectra = 0;
		int number_glyco_spectra = 0;
		
		if (print_output || draw_figures) {
		   cout << "key: \n";
		   cout << "  " << circle << " : Hexose\n";
		   cout << "  " << square << " : HexNac\n";
		   cout << "  " << triangle << " : dHex\n";
		   cout << "  " << diamond << " : NeuAc\n";
		   cout << "  X : 0,2X0 cleavage\n\n";
		}
		
		// open file for results storage
		stringstream error_tol;
		stringstream peak_intensity_tol;
		error_tol << error_tolerance;
		peak_intensity_tol << peak_intensity_tolerance;
		
		cout << "Processing... " << mgf_fname << "\n\n";
		MGF my_mgf(mgf_file, peak_intensity_tolerance);
		
		spectrum initialise;
		
		int count_spectra = 0;
		for (auto& this_spectrum : *my_mgf.get_spectra()) {
			struct spect glyco_spectrum;
			initialise.copy(&glyco_spectrum, this_spectrum);
			
			cout << "\tspectrum " << this_spectrum->title << "\n";
			//process for min mz
			initialise.remove_peaks_below_mz(&glyco_spectrum, min_mz_for_glycan);
			
			GlycanDeNovo gn = GlycanDeNovo(&glyco_spectrum, error_tolerance);
			PeptideDeNovo pn = PeptideDeNovo(this_spectrum, error_tolerance);
			
			if (pn.best_size() > min_peptide_length) {//(gn.non_isoshift_size() > min_glycan_length) && pn.best_size() > min_peptide_length) {
				cout << "\nMatches:\n\n";
				number_glyco_spectra += 1;
				
				if (draw_figures) {
					draw_annotated_spectrum(this_spectrum, pn.dn->best_graph, gn.dn->best_graph);
					draw_glycan_graph(this_spectrum, gn.dn->merged_isos_best_graph);
					remove("gnuplot_instructions.gpi");
					remove("gnuplot_instructions_glycan.gpi");
					remove("gnuplot_instructions_peptide.gpi");
					cout << "m/z range: " << gn.dn->start_node() << gn.dn->end_node() << "\n\n";
				}
				
				if (print_output) {
					
					cout << "Path consistent with peptide sequence (or reverse sequence): \n";
					pn.display();
					cout << "\nTree consistent with glycan graph: \n";
					gn.display();
					
					cout << "precursor charge = " << this_spectrum->charge << "\n";
					cout << "precursor m/z and intensity = " << this_spectrum->pepMass << "\n";
					cout << "RT in min = " << (this_spectrum->rtInSeconds / 60.0) << "\n";
					cout << "scan #= " << this_spectrum->scans << "\n";
					cout << "glycan fragmentation series charge = " << gn.dn->best_charge_state << "\n";
					
					cout << "Would you like to save this to a file?\n";
					string choice;
					while (true){
						cin >> choice;
						if (choice == "y") {
							string outfile_name = mgf_fname.substr(0, mgf_fname.length() - 4) + "_" + error_tol.str() + "Da" + peak_intensity_tol.str() + "int" + to_string(min_glycan_length) + "g.txt";
							ofstream outfile (outfile_name);
							outfile << fixed << setprecision(print_precision);
							outfile << this_spectrum->charge << " " << this_spectrum->pepMass << this_spectrum->rtInSeconds / 60.0 << " " << gn.dn->best_charge_state << " " << gn.dn->start_node() << " " << gn.dn->end_node() << "\n";
							// if output file is empty (no results saved), then remove the file from system storage
							ifstream check_outfile (outfile_name);
							if (check_outfile.peek() == ifstream::traits_type::eof())
								remove(outfile_name.c_str());
							else
								outfile.close();
							break;
						} else if (choice == "n")
							break;
						else {
							cout << "Check ya fingers and input y/n \n";
							choice = "";
						}
					}
					
					if (count_spectra < my_mgf.numSpectra() - 1)
						cout << "\ncontinued searching... \n\n";
					else
						cout << "\n* Finished searching all spectra *\n\n";
					++count_spectra;
				}
			}
			
			delete gn.dn->merged_isos_best_graph;
			delete gn.dn->best_graph;
			delete pn.dn->merged_isos_best_graph;
			delete pn.dn->best_graph;
			
		}
		for (auto& this_spectrum : *my_mgf.get_spectra())
			delete this_spectrum;
		
		mgf_file.close();
		
	} else
		cout << "\n   Incorrect number of arguments; 8 required, you gave " << argc << "\n\n";

}
