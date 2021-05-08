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

#ifndef DENOVO_H
#define DENOVO_H
#define NUM_CHARGES 3

#include<list>

struct path_pair{
	double node;
	int length;
};
struct edge{
	double u;
	double v;
	char attr;
};

#include "DiGraph.h"

void capture_drawing_data(struct spect *spectrum, Dict<double,double> *edge_series, DiGraph *best_graph) {
	for (auto& node : best_graph->Graph.nodes) {
		for (auto& _succ : best_graph->Graph.succ.get(node.first)->second) {
			edge_series->insert(node.first,spectrum->intensity[node.first]);
			edge_series->insert(_succ.first,spectrum->intensity[_succ.first]);
		}
	}
}

void write_drawing_data(string gpi, map<double,double> *intensity, Dict<double,double> *edge_series) {
	ofstream gnuplot_instructions;
	gnuplot_instructions.open(gpi);
	if (intensity != NULL)
		for (auto& mz : *intensity)
			gnuplot_instructions << mz.first << " " << mz.second << "\n";
	else
		for (auto& mz : *edge_series)
			gnuplot_instructions << mz.first << " " << mz.second << "\n";
	gnuplot_instructions.close();
}

void draw_annotated_spectrum(struct spect *spectrum, DiGraph *best_peptide, DiGraph *best_glycan) {
	string gpi = "gnuplot_instructions.gpi";											// plot spectrum
	write_drawing_data(gpi, &spectrum->intensity, NULL);
	
	string gpig = "gnuplot_instructions_glycan.gpi";
	Dict<double,double> glycan_edge_series;
	capture_drawing_data(spectrum, &glycan_edge_series, best_glycan);					// record glycan values
	write_drawing_data(gpig, NULL, &glycan_edge_series);
	
	string gpip = "gnuplot_instructions_peptide.gpi";
	Dict<double,double> peptide_edge_series;	
	capture_drawing_data(spectrum, &peptide_edge_series, best_peptide);
	write_drawing_data(gpip, NULL, &peptide_edge_series);
	
	FILE *plotHandle = popen("gnuplot -persist", "w");
	assert(plotHandle != NULL);
	fprintf(plotHandle," set terminal png; set key off; set output \"Spectrum_%s.png\"; set tics scale -0.75; set xtics 500; plot \"gnuplot_instructions.gpi\" with impulses lw 1 lc \"black\", \\\n", &spectrum->title[0]);
	fprintf(plotHandle,"\"gnuplot_instructions_glycan.gpi\" with impulses lw 3 lc \"red\", \\\n");
	fprintf(plotHandle,"\"gnuplot_instructions_peptide.gpi\" with impulses lw 3 lc \"blue\"\n");
	
	fflush(plotHandle);
	pclose(plotHandle);
	
}

void draw_glycan_graph(struct spect *spectrum, DiGraph *best_glycan) {
	ofstream gnuplot_instructions;
	gnuplot_instructions.open("gnuplot_instructions.gpi");
	
	gnuplot_instructions << fixed << setprecision(7);
	
	gnuplot_instructions << "title " << spectrum->title << "\n";
	gnuplot_instructions << "start " << best_glycan->start_node() << "\n";
	
	for (auto& node : best_glycan->Graph.nodes)
		gnuplot_instructions << node.first << "\n";
	
	for (auto& node : best_glycan->Graph.nodes) {
		for (auto& succ : best_glycan->Graph.succ.get(node.first)->second)
			gnuplot_instructions << node.first << " " << succ.first << " " << succ.second << "\n";
	}
	gnuplot_instructions.close();
	
	system("python2 plot_glycan.py");
}

class DeNovo {
	
	struct spect *current_spectrum;
	double max_distance;
	struct edge *deltas;
	double possible_charges[NUM_CHARGES] = {1.0, 2.0, 3.0};
	int mzSortedSize;
	int num_deltas;
	bool path;
	
	public:
		DiGraph *merged_isos_best_graph;
		DiGraph *best_graph;
		int best_charge_state;
		double _start_node, _end_node;
		
		typedef struct spect default_spect;
		
		DeNovo(struct spect *_current_spectrum, float _max_distance = 0.0, int num_of_deltas = 0, char _attr[] = {}, double _deltas[] = {}, bool _path = false, bool init = false){
			num_deltas = num_of_deltas*NUM_CHARGES;
			struct edge delta_array[num_deltas];
			deltas = &delta_array[0];
			max_distance = _max_distance;
			if (init) {
				current_spectrum = _current_spectrum;
				int iso_count[NUM_CHARGES] = {0, -1, 1};
				path = _path;
				
				for (int i = 0, j = 0; i < num_of_deltas; ++i) {
					for (int iso_count_val = 0; iso_count_val < NUM_CHARGES; ++iso_count_val, ++j) {
						double a = _deltas[i] + iso_count[iso_count_val]*1.003;
						struct edge new_delta;
						new_delta.u = a;
						new_delta.attr = _attr[i];
						deltas[j] = new_delta;
					}
				}
				init_best_graph_and_charge_state();
			}
		}
		
		~DeNovo() {}
		
		void add_isoshift_edges(DiGraph* dg) {
			list<edge> edges_to_add;
			// for all sibling nodes B, C where n --> B and n --> C, add
			// edges from:
			// 		-B to all successors of C
			// 		-C to all successors of B
			// 		-all predecessors of B to C
			// 		-all predecessors of C to B
			// this ensures the nodes are treated properly (you can go in one and come out the other)
			
			for (auto& a : dg->Graph.nodes) {
				Dict<double,char>* successors = &dg->Graph.succ.get(a.first)->second;
				for (auto& succB : *successors) {
					for (auto& succC : *successors) {
						
						if (succB.second == succC.second) {
							if (succB.first == succC.first)
								continue;
							for (auto& desc : dg->Graph.succ.get(succC.first)->second) {
								struct edge new_edge;
								new_edge.u = succB.first;
								new_edge.v = desc.first;
								new_edge.attr = desc.second;
								edges_to_add.push_back(new_edge);
							}
							
							for (auto& anc : dg->Graph.pred.get(succC.first)->second) {
								struct edge new_edge;
								new_edge.u = anc.first;
								new_edge.v = succB.first;
								new_edge.attr = dg->Graph.succ.get(anc.first)->second.get(succC.first)->second;
								edges_to_add.push_back(new_edge);
							}
						}
					}
				}
			}
			for (auto& a : dg->Graph.nodes) {
				Dict<double,char>* predecessors = &dg->Graph.pred.get(a.first)->second;
				for (auto& predB : *predecessors)
					for (auto& predC : *predecessors) {
						
						if (dg->Graph.succ.get(predB.first)->second.get(a.first)->second == dg->Graph.succ.get(predC.first)->second.get(a.first)->second) {
							if (predB.first == predC.first)
								continue;
							
							for (auto& desc : dg->Graph.succ.get(predC.first)->second) {
								struct edge new_edge;
								new_edge.u = predB.first;
								new_edge.v = desc.first;
								new_edge.attr = desc.second;
								edges_to_add.push_back(new_edge);
							}
							for (auto& anc : dg->Graph.pred.get(predC.first)->second) {
								struct edge new_edge;
								new_edge.u = anc.first;
								new_edge.v = predB.first;
								new_edge.attr = dg->Graph.succ.get(anc.first)->second.get(predC.first)->second;
								edges_to_add.push_back(new_edge);
							}
						}
					}
			}
			for (auto& each_edge : edges_to_add)
				dg->add_edge(each_edge.u, each_edge.v,	each_edge.attr);
		}
		
		double peak_nearest_mz_helper(double mz, int low, int high, double mzSorted[]) {
			// note: could be adapted to use a hash table and windows for
			// greater speed (particularly when error tolerance is low);
			// also, all peaks within the error tolerance could be found,
			// not simply the closest.
			if (low == high)
				return mzSorted[low];
			if (low + 1 == high) {
				// arbitrarily resolves ties
				if (abs(mzSorted[low] - mz) < abs(mzSorted[high] - mz)) {
					return mzSorted[low];
				} else {
					return mzSorted[high];
				}
			}
			int mid = (int) (low + (high - low + 1) / 2);
			if (mz < mzSorted[mid])
				high = mid;
			else if (mz > mzSorted[mid])
				low = mid;
			else{
				return mz;
			}
			return peak_nearest_mz_helper(mz, low, high, mzSorted);
		}
		
		double peak_nearest_mz(double mz, double mzSorted[]) {
			// mzSorted is composed of sorted, unique elements
			assert(mzSortedSize > 0);
			if (mzSortedSize == 1)
				return mzSorted[0];
			int low = 0, high = mzSortedSize - 1;
			return peak_nearest_mz_helper(mz, low, high, mzSorted);
		}
		
		void build_graph(DiGraph* dg, double _charge) {
			mzSortedSize = current_spectrum->intensity.size();
			double mzSorted[mzSortedSize];
			int i = 0;
			for (map<double,double>::iterator intensity_vals = current_spectrum->intensity.begin(); intensity_vals != current_spectrum->intensity.end(); ++intensity_vals, ++i) {
				mzSorted[i] = intensity_vals->first;
			}
			for(int mz = 0; mz < mzSortedSize; ++mz) {
				for(int i = 0; i < num_deltas; ++i) {
					double change = deltas[i].u / _charge;
					double neighborMz = peak_nearest_mz(mzSorted[mz] + change, mzSorted);
					double distance = abs(mzSorted[mz] + change - neighborMz);
					if (distance < max_distance && neighborMz > mzSorted[mz]) {
						dg->add_edge(mzSorted[mz], neighborMz, deltas[i].attr);
					}
				}
			}
			add_isoshift_edges(dg);
		}
		
		void init_best_graph_and_charge_state() {
			list<DiGraph*> results_for_all_charges;
			for(int charge = 0; charge < sizeof(possible_charges)/sizeof(possible_charges[0]); ++charge) {
				DiGraph dg;
				try {
					build_graph(&dg, possible_charges[charge]);
					if (dg.size() == 0) {
						results_for_all_charges.push_back(new DiGraph);
						throw "";
					}

					DiGraph *best_graph_for_charge = new DiGraph;
					if (!path) {
						maximal_directed_connected_graph(best_graph_for_charge, &dg);
					} else {
						maximal_paths_graph(best_graph_for_charge, &dg);
					}
					best_graph_for_charge->Graph.charge = possible_charges[charge];
					
					int _best_graph = best_graph_for_charge->Graph.nodes.size();
					if (_best_graph > 0)
						results_for_all_charges.push_back(best_graph_for_charge);
				} catch (char const *e) {}
			}
			if (results_for_all_charges.size() > 0) {
				best_graph_all_charges(results_for_all_charges);
				best_charge_state = best_graph->Graph.charge;
				
				for (auto& each_result : results_for_all_charges) {
					if (each_result != best_graph)
						delete each_result;
				}
				merged_isos_best_graph = new DiGraph;
				merge_isos();
				_start_node = merged_isos_best_graph->start_node();
				_end_node = merged_isos_best_graph->end_node();
			}
		}
		
		void maximal_directed_connected_graph(DiGraph *best_undirected_graph_sub, DiGraph* dg) {
			list<DiGraph*> connected_graphs;
			dg->connected_component_subgraphs(&connected_graphs,dg);
			
			//maximal subgraph by number of nodes
			DiGraph *best_undirected_graph = dg->max_graph(&connected_graphs);
			//best_undirected_graph->print_graph();

			//now subgraph of dg
			list<double> best_undirected_graph_nodes;
			for (auto& _node : best_undirected_graph->Graph.nodes)
				best_undirected_graph_nodes.push_back(_node.first);
			
			dg->subgraph_directed(best_undirected_graph_sub, &best_undirected_graph->Graph.nodes);
			
			for (auto& each_graph : connected_graphs) {
				//each_graph->print_graph();
				delete each_graph;
			}
		}
			
		void maximal_paths_graph(DiGraph *best_undirected_graph, DiGraph* dg) {
			list<DiGraph*> connected_paths;
			dg->longest_paths(&connected_paths);
			best_undirected_graph->deep_copy(dg->max_graph(&connected_paths));
			for (auto& each_graph : connected_paths)
				delete each_graph;
		}
		
		void best_graph_all_charges(list<DiGraph*> results_for_all_charges){
			int max_size = 0;
			for (list<DiGraph*>::iterator _graph = results_for_all_charges.begin(); _graph != results_for_all_charges.end(); _graph++) {
				int graph_size = (*_graph)->Graph.nodes.size();
				if (max_size <= graph_size) {
					max_size = graph_size;
					best_graph = (*_graph);
				}
			}
		}

		void merge_isos() {
			Dict<double,double> merged_nodes;
			
			for (auto& a : best_graph->Graph.nodes) {
				for (auto& b : best_graph->Graph.nodes) {
					if (a.first == b.first)
						continue;
					if (has_merge_nodes(best_graph, &merged_nodes, a.first, b.first)) {
						continue;
					}
				}
			}
			
			merged_isos_best_graph->deep_copy_by_nodes(best_graph);
			for (auto& a : merged_nodes) {
				if (!merged_isos_best_graph->Graph.nodes.contains(a.first) || !merged_isos_best_graph->Graph.nodes.contains(a.second))
					continue;
				for (auto &pred : merged_isos_best_graph->Graph.pred.get(a.second)->second) {
					dict_item<double,Dict<double,char>> *predecessor = merged_isos_best_graph->Graph.succ.get(pred.first);
					dict_item<double,char> *_predecessor = predecessor->second.get(a.second);
					if (predecessor != NULL && _predecessor != NULL)
						merged_isos_best_graph->add_edge(pred.first, a.first, _predecessor->second);
				}
				for (auto &succ: merged_isos_best_graph->Graph.succ.get(a.second)->second) {
					dict_item<double,Dict<double,char>> *successor = merged_isos_best_graph->Graph.succ.get(succ.first);
					dict_item<double,char> *_successor = successor->second.get(a.second);
					if (successor != NULL && _successor != NULL)
						merged_isos_best_graph->add_edge(a.first, succ.first, _successor->second);
				}
				merged_isos_best_graph->remove_node(a.second);
			}
			//merged_isos_best_graph->print_graph();
		}
		
		bool has_merge_nodes(DiGraph* dg, Dict<double,double>* merged_nodes, double a, double b) {
			for (auto& p : dg->Graph.pred.get(a)->second)
				for (auto& q : dg->Graph.pred.get(b)->second)
					if (p.first == q.first || merged_nodes->contains_pair(p.first, q.first)) {
						if (dg->Graph.succ.get(p.first)->second.get(a)->second == dg->Graph.succ.get(q.first)->second.get(b)->second) {
							merged_nodes->insert(a,b,true);
							merged_nodes->insert(b,a,true);
							return true;
						}
					}
			for (auto& p : dg->Graph.succ.get(a)->second)
				for (auto& q : dg->Graph.succ.get(b)->second)
					if (p.first == q.first || merged_nodes->contains_pair(p.first, q.first)) {
						if (dg->Graph.succ.get(a)->second.get(p.first)->second == dg->Graph.succ.get(b)->second.get(q.first)->second) {
							merged_nodes->insert(a,b,true);
							merged_nodes->insert(b,a,true);
							return true;
						}
					}
			return false;
		}
		
		int non_isoshift_size() {
			if (merged_isos_best_graph == NULL)
				return 0;
			else
				return merged_isos_best_graph->Graph.nodes.size();
		}
		
		int best_graph_size() {
			if (best_graph == NULL)
				return 0;
			else
				return best_graph->Graph.nodes.size();
		}
		
		double start_node() {
			return _start_node;
		}
		
		double end_node() {
			return _end_node;
		}
		
};

class GlycanDeNovo {
	public:
		DeNovo *dn;
		double iso_shift_constant = 1.003;
		GlycanDeNovo(struct spect *current_spectrum, float max_distance) {
			char attr[] = {circle, square, triangle, diamond};
			double deltas[] = {162.0528, 203.0794, 146.0579, 291.0954};
			
			dn = new DeNovo(current_spectrum, max_distance, sizeof(deltas)/sizeof(deltas[0]), attr, deltas, false, true);
		}
		
		~GlycanDeNovo() {
			delete dn;
		}
		
		int non_isoshift_size() {
			return dn->non_isoshift_size();
		}
		
		void display() {
			cout << fixed << setprecision(print_precision);
			cout << "charge state: " << dn->best_charge_state << "\n";
			cout << "graph (with merged isotope peaks) is in m/z range: " << dn->_start_node << " " << dn->_end_node << "\n";
			display_graph();
		}
		
		void display_graph() {
			set<double> visited_nodes;
			// can also print the graph where redundant paths from isotopic
			// shifts are merged by commenting out the second line and
			// commenting the first:
			DiGraph *best_graph = dn->best_graph;
			// DiGraph best_graph = dn->merged_isos_best_graph;
			for (auto& n : best_graph->Graph.nodes) {
				if (best_graph->Graph.pred.get(n.first)->second.size() == 0)
					display_graph_helper(n.first, &visited_nodes);
			}
			cout << "\n";
		}
		
		void display_graph_helper(double root, set<double> *visited_nodes, int depth = 0) {
			cout << root << "\n";
			(*visited_nodes).insert(root);
			for (auto& n : dn->best_graph->Graph.succ.get(root)->second)
				if (*(*visited_nodes).find(n.first) != n.first) {
					indent(depth+1);
					cout << " + " << dn->best_graph->Graph.succ.get(root)->second.get(n.first)->second << " --> ";
					display_graph_helper(n.first, visited_nodes, depth+1);
				} else {
					indent(depth+1);
					char attr = dn->best_graph->Graph.succ.get(root)->second.get(n.first)->second;
					if (attr == 'I')
						cout << " + " << "I/L" << " " << n.first << "\n";
					else cout << " + " << attr << " " << n.first << "\n";
				}
		}
		
		void indent(int depth) {
			for (int i = 0; i < depth; ++i)
				cout << "     ";
		}
		
};

class PeptideDeNovo {
	public:
		DeNovo *dn;
		PeptideDeNovo(struct spect *current_spectrum, double max_distance) {
			char attr[] = {'A','R','N','D','C','E','Q','G','H','I','K','M','F','P','S','T','W','Y','V'};
			double deltas[] = {71.03711, 156.10111, 114.04293, 115.02694, 103.00919, 129.04259, 128.05858, 57.0519, 137.05891, 113.08406, 128.09496, 131.04049, 147.06841, 97.05276, 87.03203, 101.04768, 186.07931, 163.06333, 99.06841};
			
			dn = new DeNovo(current_spectrum, max_distance, sizeof(deltas)/sizeof(deltas[0]), attr, deltas, true, true);
		}
		
		~PeptideDeNovo() {
			delete dn;
		}
		
		int best_size() {
			return dn->best_graph_size();
		}
		
		void display() {
			cout << "charge state: " << dn->best_charge_state << "\n";
			display_graph();
		}
		
		void display_graph() {
			map<double, double> ordered_nodes;
			for (auto& _node : dn->best_graph->Graph.nodes)
				ordered_nodes[_node.first] = _node.first;
			for (auto& in_order : ordered_nodes)
				if (dn->best_graph->Graph.succ.get(in_order.first)->second.size() > 0) {
					char attr = dn->best_graph->Graph.succ.get(in_order.first)->second.begin()->second;
					if (attr == 'I')
						cout << in_order.first << " " << dn->best_graph->Graph.succ.get(in_order.first)->second.begin()->first << " " << "I/L" << "\n";
					else
						cout << in_order.first << " " << dn->best_graph->Graph.succ.get(in_order.first)->second.begin()->first << " " << attr << "\n";
				}
		}
		
};

#endif
