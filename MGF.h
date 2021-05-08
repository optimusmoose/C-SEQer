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

#ifndef MGF_H
#define MGF_H

#include<cstring>

struct spect{
	string title;
	string scans;
	int charge;
	string pepMass;
	double rtInSeconds;
	map<double, double> intensity;
};

class spectrum{
	public:
		void read_spectra(ifstream& file, struct spect *newSpectrum, double _percent_peak_threshold){
			for (string line; getline(file, line);) {
				int lineLen = line.length();
				if (lineLen == 0)
					continue;
				if (line.compare(0, 1, "#") == 0)
					continue;
				if (line.compare(0, 10, "BEGIN IONS") == 0)
					continue;
				else if (line.compare(0, 8, "END IONS") == 0)
					break;
				else if (line.compare(0, 6, "TITLE=") == 0) {
					//cout << line.substr(6, lineLen - 1) << "\n";
					newSpectrum->title = line.substr(6, lineLen - 1);
				} else if (line.compare(0, 6, "SCANS=") == 0) {
					//cout << line.substr(6, lineLen - 1) << "\n";
					newSpectrum->scans = line.substr(6, lineLen - 1);
				} else if (line.compare(0, 7, "CHARGE=") == 0) {
					//cout << line.substr(7, lineLen - 1) << "\n";
					newSpectrum->charge = stoi(line.substr(7, lineLen - 1));
				} else if (line.compare(0, 8, "PEPMASS=") == 0) {
					//cout << line.substr(8, lineLen - 1) << "\n";
					newSpectrum->pepMass = line.substr(8, lineLen - 1);
					newSpectrum->pepMass[newSpectrum->pepMass.length()-1] = ' ';
				} else if (line.compare(0, 7, "RTINSEC") == 0) {
					//cout << line.substr(12, lineLen - 1) << "\n";
					newSpectrum->rtInSeconds = stof(line.substr(12, lineLen - 1));
				} else {
					int split = line.find(" ");
					if (split == -1)
						split = line.find("\t");
					if (split == -1)
						continue;
					string mz_str = line.substr(0, split), intensity_val_str = line.substr(split + 1, lineLen - 1);
					double mz = stod(mz_str);
					double intensity_val = stod(intensity_val_str);
					//set precision for printing
					//cout << fixed;
					//cout << setprecision(print_precision);
					//cout << mz << " " << intensity_val << "\n";
					newSpectrum->intensity.insert( pair<double, double>(mz, intensity_val) );
					//newSpectrum->intensity[mz] = intensity_val;
				}
			}
			init_spectrum(newSpectrum, _percent_peak_threshold);
		}
		
		void init_spectrum(struct spect *_spectToInit, double percent_peak_threshold) {
			double max_intensity = 0.0;
			for (map<double,double>::iterator intensity_vals = _spectToInit->intensity.begin(); intensity_vals != _spectToInit->intensity.end(); ++intensity_vals)
				if (intensity_vals->second > max_intensity)
					max_intensity = intensity_vals->second;
			// copy matching entries to new map
			//cout << fixed;
			//cout << setprecision(print_precision);
			//cout << "\n\t\t" << percent_peak_threshold << ", " << max_intensity << ", " << max_intensity / (1/percent_peak_threshold) << "\n";
			double _percent_peak_threshold = (1/percent_peak_threshold);
			map<double,double> newMap;
			for (map<double,double>::iterator intensity_vals = _spectToInit->intensity.begin(); intensity_vals != _spectToInit->intensity.end(); ++intensity_vals)
				if (intensity_vals->second > max_intensity / _percent_peak_threshold)
					newMap[intensity_vals->first] = intensity_vals->second;
			_spectToInit->intensity = newMap;
		}
		
		void remove_peaks_below_mz(struct spect *_spectToInit, double mz_threshold) {
			map<double,double> newMap;
			for (map<double,double>::iterator intensity_vals = _spectToInit->intensity.begin(); intensity_vals != _spectToInit->intensity.end(); ++intensity_vals)
				if (intensity_vals->first > mz_threshold)
					newMap[(intensity_vals)->first] = (intensity_vals)->second;
			_spectToInit->intensity = newMap;
		}
	
		void copy(struct spect *copyInto, struct spect *copyFrom) {
			// only needs to copy the intensity values
			for (map<double,double>::iterator intensity_vals = copyFrom->intensity.begin(); intensity_vals != copyFrom->intensity.end(); ++intensity_vals)
				copyInto->intensity.insert( pair<double, double>(intensity_vals->first, intensity_vals->second) );
		}
		
		void print_spectrum(struct spect *currentSpectrum) {
			cout << fixed << setprecision(print_precision);
			if (currentSpectrum->title != "")
				cout << "\n\tTitle: " << currentSpectrum->title;
			if (currentSpectrum->scans != "")
				cout << "\n\tscans: " << currentSpectrum->scans;
			if (currentSpectrum->charge != 0)
				cout << "\n\tCharge: " << currentSpectrum->charge;
			if (currentSpectrum->pepMass != "")			
				cout << "\n\tPep Mass: " << currentSpectrum->pepMass;
			if (currentSpectrum->rtInSeconds > 0.0)
				cout << "\n\trt in seconds: " << currentSpectrum->rtInSeconds;
			cout << "\n";
			for (map<double,double>::iterator intensity_vals = currentSpectrum->intensity.begin(); intensity_vals != currentSpectrum->intensity.end(); ++intensity_vals)
				cout << intensity_vals->first << "\t\t" << intensity_vals->second << '\n';
			cout << "\t*** end spectrum ***\n";
		}
};

class MGF {
	// instance variables
	int num_spectra = 0;
	vector<struct spect*> all_spectra;
	spectrum new_spectrum;

	public:
		spectrum test;
		// constructor
		MGF(ifstream& file, double percent_peak_threshold){
			int lineCounter = 0;
			int spectrumCounter = 0;
			while (!file.eof()) {
				struct spect *temp_spectrum = new (struct spect);
				new_spectrum.read_spectra(file, temp_spectrum, percent_peak_threshold);
				if (temp_spectrum->intensity.size() > 0) {
					all_spectra.emplace_back(temp_spectrum);
					++num_spectra;
				} else {
					delete temp_spectrum;
				}
			}
		}	
		
		~MGF() {}
		
		// display spectra to check for correctness
		void print_mgf_spectra() {
			for(int i = 0; i < all_spectra.size(); i++) {
				cout << fixed;
				cout << setprecision(print_precision);
				if (all_spectra[i]->title != "")
					cout << "\n\tTitle: " << all_spectra[i]->title;
				if (all_spectra[i]->scans != "")
					cout << "\n\tscans: " << all_spectra[i]->scans;
				if (all_spectra[i]->charge != 0)
					cout << "\n\tCharge: " << all_spectra[i]->charge;
				if (all_spectra[i]->pepMass != "")
					cout << "\n\tPep Mass: " << all_spectra[i]->pepMass;
				if (all_spectra[i]->rtInSeconds > 0.0)
					cout << "\n\trt in seconds: " << all_spectra[i]->rtInSeconds;
				cout << "\n";
				for (map<double,double>::iterator intensity_vals = all_spectra[i]->intensity.begin(); intensity_vals != all_spectra[i]->intensity.end(); ++intensity_vals)
					cout << intensity_vals->first << "\t\t" << intensity_vals->second << '\n';
				cout << "\t*** end spectrum ***\n";
			}
		}
		
		vector<struct spect*>* get_spectra(){
			return &all_spectra;
		}
		
		int numSpectra(){
			return num_spectra;
		}
};

#endif