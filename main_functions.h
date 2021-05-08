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

#ifndef MAIN_FUNCTIONS
#define MAIN_FUNCTIONS

const char circle = 'H';
const char square = 'N';
const char diamond = 'S';
const char triangle = 'F';

bool boolean(string parameter) {
	if (parameter.compare("0") == 0)
		return false;
	else if (parameter.compare("1") == 0)
		return true;
	else 
		throw "Parameters 6 and 7 must be 0 or 1";
}

void print_parameters(double error_tolerance, double peak_intensity_tolerance, int min_peptide_length, int min_glycan_length, double min_mz_for_glycan, bool draw_figures) {
	cout << "\nYou called SweetSEQer with these parameters:\n";
		cout << "\tepsilon = " << error_tolerance << "\n";
		cout << "\tlambda = " << peak_intensity_tolerance << "\n";
		cout << "\tp = " << min_peptide_length << "\n";
		cout << "\tg = " << min_glycan_length << "\n";
		cout << "\tSearching for glycan fragments at or above " << min_mz_for_glycan << " m/z\n";
		if (draw_figures == true)
			cout << "\tDrawing figures for spectra of note\n\n";
}

void print_usage(char const *exception) {
	cout << "\nERROR: " << exception;
	cout << "\n\n usage: ./SweetSEQer.py <epsilon> <lambda> <p> <g> <glycan_min_m/z> <plot_results> <MGF_file_1> [MGF_file_2...]";
	cout << "\n\tepsilon: maximum absolute error between predicted and actual peak (in m/z)";
	cout << "\n\tlambda: minimum peak intensity (relative to maximum peak intensity)";
	cout << "\n\tp: minimum length of inferred peptide sequence";
	cout << "\n\tg: minimum size of inferred glycan graph (including isotope peaks, not shown in graph)";
	cout << "\n\ttau: minimum m/z value to include in glycan";
	cout << "\n\tprint_text_output: either 0 (doesn\'t print text results) or 1 (prints text results for matches)";
	cout << "\n\tplot_results: either 0 (doesn\'t plot results) or 1 (interactively plots each matching spectrum)";
	cout << "\n\tMGF_file_1 [MGF_file_2...]: paths to an MGF files to process \n\n";
}

void capture_parameters(char *args[], double* error_tolerance, double* peak_intensity_tolerance, int* min_peptide_length, int* min_glycan_length, double* min_mz_for_glycan, bool* print_output, bool* draw_figures, string* mgf_fname, ifstream& mgf_file) {
	try {																						// assignment of parameters and exception handling for acceptable values
		*error_tolerance = stof(args[1]);
		if (*error_tolerance <= 0)
			throw "epsilon must be > 0";
		
		*peak_intensity_tolerance = stof(args[2]);
		if (*peak_intensity_tolerance <= 0)
			throw "lambda must be > 0";
		
		*min_peptide_length = stoi(args[3]);
		if (*min_peptide_length < 1)
			throw "p must be >= 1";
		
		*min_glycan_length = stoi(args[4]);
		if (*min_glycan_length < 1)
			throw "g must be >= 1";
		
		*min_mz_for_glycan = stof(args[5]);
		if (*min_mz_for_glycan < 0)
			throw "tau must be >= 0";
		
	} catch (char const* e) {
		print_parameters(*error_tolerance, *peak_intensity_tolerance, *min_peptide_length, *min_glycan_length, *min_mz_for_glycan, false);
		print_usage(e);
		exit (EXIT_FAILURE);
	}
	
	try {
		*print_output = boolean(args[6]);
		*draw_figures = boolean(args[7]);
	} catch (char const* e) {
		print_parameters(*error_tolerance, *peak_intensity_tolerance, *min_peptide_length, *min_glycan_length, *min_mz_for_glycan, *draw_figures);
		print_usage(e);
		exit (EXIT_FAILURE);
	}
	
	try {
		*mgf_fname = args[8];
		int len = (*mgf_fname).length();
		if ((*mgf_fname).compare(len - 4, len, ".mgf") != 0)
			throw "file name must end with '.mgf'";
		
		mgf_file.open(*mgf_fname);
		if (mgf_file.is_open() != 1)
			throw "specified mgf file could not be opened";
		
	} catch (char const* e) {
		print_usage(e);
		exit (EXIT_FAILURE);
	}
}

#endif