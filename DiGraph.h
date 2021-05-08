/*    
	Copyright (C) 2021-2035
	This minimal library was written in C++ based on the original python library for NetworkX as written by:
		Aric Hagberg <hagberg@lanl.gov>
		Dan Schult <dschult@colgate.edu>
		Pieter Swart <swart@lanl.gov>
    All rights reserved.
    
	License: BSD (Berkeley Software Distribution)
	
	Author for C++: Chris Burgoyne <christopher.burgoyne@umontana.edu>
	***Base class for directed graphs for the SweetSEQer algorithm written in C++***
*/

#ifndef DiGraph_H
#define DiGraph_H

#include<set>
#include "Dictionary.h"	// custom library for Python(2.7.6) dictObject

struct graph{
	Dict<double, double> nodes;				// nodes
	Dict<double, Dict<double, char>> pred;	// predecessor
	Dict<double, Dict<double, char>> succ;	// successor
	int charge;
	int longest_path_length;
};

class DiGraph {
	public:
	struct graph Graph;
	
	DiGraph(){}
	
	~DiGraph(){}
	
    void add_edge(double u, double v, char attr) {
		
		if (!Graph.nodes.contains(u)) {
			Graph.nodes.insert(u,u);
			Graph.pred.insert(u);
			Graph.succ.insert(u);
		}
		if (!Graph.nodes.contains(v)) {
			Graph.nodes.insert(v,v);
			Graph.pred.insert(v);
			Graph.succ.insert(v);
		}
		Graph.succ.get(u)->second.insert(v, attr);
		Graph.pred.get(v)->second.insert(u, attr);
	}
	
	void remove_node(double node) {
		Graph.nodes.erase(node);
		Graph.pred.erase(node);
		Graph.succ.erase(node);
		
		for (Dict<double,Dict<double,char>>::Iterator _succ = Graph.succ.begin(); _succ != Graph.succ.end(); ++_succ)
			if (_succ->second.contains(node)) {
				//cout << "\n***" << _succ->first << " " << node << " " << _succ->second.contains(node) << " ***\n";
				_succ->second.erase(node);
			}
		for (Dict<double,Dict<double,char>>::Iterator _pred = Graph.pred.begin(); _pred != Graph.pred.end(); ++_pred)
			if (_pred->second.contains(node)) _pred->second.erase(node);
	}
	
	bool has_edge(double u, double v, char attr) {
		if (Graph.succ.contains(u)) { 
			Dict<double,char> *_succ = &Graph.succ.get(u)->second;
			if (_succ->contains(v))
				if (_succ->get(v)->second == attr)
					return true;
		}
		return false;
	}
	
	static bool compare_conn_components(const list<double> first, const list<double> second) {
		if (first.size() > second.size())
			return true;
		return false;
	}
	
	void connected_component_subgraphs(list<DiGraph*>* ccs, DiGraph *G) {
		DiGraph undirected;
		to_undirected(&undirected, G);
		
		list<list<double>> cc;
		connected_components(&cc, &undirected);
		cc.sort(compare_conn_components);
		for (auto& c : cc)
			ccs->push_back(undirected.subgraph(&c));
	}
	
	void connected_components(list<list<double>> *components, DiGraph *G) {
		Dict<double,int> seen;
		for (auto& v : G->Graph.nodes)
			if (!seen.contains(v.first)) {
				Dict<double,int> c;
				single_source_shortest_path_length(&c, G, v.first);
				list<double> c_c;
				for (auto& _node : c) {
					c_c.push_back(_node.first);
					//cout << _node.first << " ";
				}
				//cout << "\n\n";
				components->push_front(c_c);
				seen.update(&c);
			}
		components->sort(compare_length);
	}
	
	static bool compare_length(const list<double> first, const list<double> second) {
		return first.size() >= second.size();
	}
	
	void single_source_shortest_path_length(Dict<double,int> *seen, DiGraph *G, double source) {
		int level = 0;
		Dict<double,char> *nextlevel = new Dict<double,char>;
		nextlevel->insert(source,1);
		while (nextlevel->size() != 0) {
			Dict<double,char> *thislevel = nextlevel;
			nextlevel = new Dict<double,char>;
			for (auto& v : *thislevel){
				if (!seen->contains(v.first)) {
					//cout << v.first << "*\n";
					seen->insert(v.first,level);
					nextlevel->update(&G->Graph.succ.get(v.first)->second);
					/* for (auto &item : G->Graph.succ.get(v.first)->second)
						cout << item.first << " ";
					cout << "\n";
					for (auto &item : *nextlevel)
						cout << item.first << " ";
					cout << "\n"; */
					//cout << v.first << "\n";
				}
			}
			//cout << " ***\n";
			++level;
			delete thislevel;
		}
		delete nextlevel;
	}
	
	DiGraph* subgraph(list<double> *nbunch) {
		DiGraph *subgraph = new DiGraph;
		for (auto& n : *nbunch) {
			subgraph->Graph.nodes.insert(n,n);
			subgraph->Graph.succ.insert(n);
			for (auto& nbr : Graph.succ.get(n)->second) {
				if (in_dict(nbr.first, &subgraph->Graph.succ)) {
					subgraph->Graph.succ.get(n)->second.insert(nbr.first,nbr.second);
					if (!subgraph->Graph.succ.contains(nbr.first))
						subgraph->Graph.succ.insert(nbr.first);
					subgraph->Graph.succ.get(nbr.first)->second.insert(n,nbr.second);
				}
			}
		}
		return subgraph;
	}
	
	void subgraph_directed(DiGraph *subgraph, Dict<double,double> *nbunch) {
		for (auto& n : *nbunch) {
			subgraph->Graph.nodes.insert(n.first,n.first);
			subgraph->Graph.succ.insert(n.first);
			subgraph->Graph.pred.insert(n.first);
		}
		for (auto& u : subgraph->Graph.succ) {
			for (auto& v : Graph.succ.get(u.first)->second) {
				if (in_dict(v.first, &subgraph->Graph.succ)) {
					subgraph->Graph.succ.get(u.first)->second.insert(v.first,v.second);
					subgraph->Graph.pred.get(v.first)->second.insert(u.first,v.second);
				}
			}
		}
	}
	
	bool in_dict(double key, Dict<double,Dict<double,char>> *H) {
		for (auto& _each : *H) {
			if (key == _each.first)
				return true;
			for (auto& _each_dbl : _each.second)
				if (key == _each_dbl.first)
					return true;
		}
		return false;
	}

	void to_undirected(DiGraph *undirected, DiGraph *dg) {
		list<edge> edges_to_add;
		// add nodes
		for (auto& _node : dg->Graph.nodes) {
			undirected->Graph.nodes.insert(_node.first, _node.second);
		}
		
		// add edges in both directions
		for (auto& _node : dg->Graph.nodes) {
			for (auto& succ : dg->Graph.succ.get(_node.first)->second) {
				struct edge new_edge;
				new_edge.u = _node.first;
				new_edge.v = succ.first;
				new_edge.attr = succ.second;
				edges_to_add.push_back(new_edge);
			}
		}
		
		for (auto& _edge : edges_to_add) {
			if (!in_dict(_edge.u, &undirected->Graph.succ))
				undirected->Graph.succ.insert(_edge.u);
			if (!in_dict(_edge.v, &undirected->Graph.succ))
				undirected->Graph.succ.insert(_edge.v);
			undirected->Graph.succ.get(_edge.u)->second.insert(_edge.v,_edge.attr);
			undirected->Graph.succ.get(_edge.v)->second.insert(_edge.u,_edge.attr);
		}
	}
	
	void longest_paths(list<DiGraph*>* results){
		Dict<double,int> node_to_len_of_longest_path_ending;
		list<struct path_pair> longest_path_len_and_node_list;
		
		for (auto& _node : Graph.nodes) {
			struct path_pair new_pair;
			new_pair.node = _node.first;
			new_pair.length = len_of_longest_path_ending(_node.second, _node.second, &node_to_len_of_longest_path_ending, &longest_path_len_and_node_list);
			longest_path_len_and_node_list.push_back(new_pair);
			
		}
		
		int best_len = 0;
		double tail_of_longest_path = 0.0;
		list <struct path_pair> equal_length_paths;
		for (auto& _pair : longest_path_len_and_node_list) {		// collect all paths of equal length
			if (_pair.length >= best_len) {
				best_len = _pair.length;
				struct path_pair new_pair;
				new_pair.node = _pair.node;
				new_pair.length = _pair.length;
				equal_length_paths.push_back(new_pair);
			}
		}
		for (auto& _pair : equal_length_paths) {					// find longest path based on max node value
			if (_pair.length < best_len)
				continue;
			else
				if (_pair.node > tail_of_longest_path)
					tail_of_longest_path = _pair.node;
		}
		
		// create best path graph for charge
		list<struct edge> edges;
		while (true) {
			Dict<double,char> *preds = &Graph.pred.get(tail_of_longest_path)->second;
			if (preds->size() == 0)
				break;
			--best_len;
			for (auto& n : *preds)
				if (node_to_len_of_longest_path_ending.get(n.first)->second == best_len) {
					struct edge new_edge;
					new_edge.u = n.first;
					new_edge.v = tail_of_longest_path;
					edges.push_front(new_edge);
					tail_of_longest_path = n.first;
					break;
				}
		}
		// create graph of best path for charge
		DiGraph *result = new DiGraph;
		for (auto& _edge : edges) {
			result->add_edge(_edge.u, _edge.v, Graph.succ.get(_edge.u)->second.get(_edge.v)->second);
		}
		results->push_back(result);
	}
	
	int len_of_longest_path_ending(double n, double node, Dict<double,int> *node_to_len_of_longest_path_ending, list<struct path_pair> *longest_path_len_and_node_list) {
		if (node_to_len_of_longest_path_ending->contains(n)) {
			return node_to_len_of_longest_path_ending->get(n)->second;
		}
		node_to_len_of_longest_path_ending->insert(n, 0);
		Dict<double,char> *preds = &Graph.pred.get(n)->second;
		for (auto& m : *preds) {
			int max;
			int a = node_to_len_of_longest_path_ending->get(n)->second; 
			int b = 1 + len_of_longest_path_ending(m.first, node, node_to_len_of_longest_path_ending, longest_path_len_and_node_list);
			if (a >= b) max = a;
			else max = b;
			node_to_len_of_longest_path_ending->insert(n, max);
		}
		return node_to_len_of_longest_path_ending->get(n)->second;
	}
	
	DiGraph* max_graph(list<DiGraph*>* connected_graphs) {												// what should this method return when there are multiple subgraphs of the same size?
		DiGraph *maximal_graph;
		int max_size = 0;
		for (list<DiGraph*>::iterator _graph = connected_graphs->begin(); _graph != connected_graphs->end(); ++_graph) {int graph_size = (*_graph)->Graph.nodes.size();
			if (max_size < graph_size) {
				max_size = graph_size;
				maximal_graph = *_graph;
			}
		}
		return maximal_graph;
	}
	
	void deep_copy(DiGraph *dg) {
		Graph.nodes.copy_nodes(&dg->Graph.nodes);
		for (auto& node : dg->Graph.nodes) {
			Graph.pred.insert(node.first);
			Graph.succ.insert(node.first);
			Graph.pred.get(node.first)->second.copy_nodes(&dg->Graph.pred.get(node.first)->second);
			Graph.succ.get(node.first)->second.copy_nodes(&dg->Graph.succ.get(node.first)->second);
			
		}
	}
	
	void deep_copy_by_nodes(DiGraph *dg) {
		for (auto& node : dg->Graph.nodes) {
			Graph.nodes.insert(node.first,node.first);
			Graph.succ.insert(node.first);
			Graph.pred.insert(node.first);
		for (auto& succ : dg->Graph.succ.get(node.first)->second)
			Graph.succ.get(node.first)->second.insert(succ.first, succ.second);
		for (auto& pred : dg->Graph.pred.get(node.first)->second)
			Graph.pred.get(node.first)->second.insert(pred.first, pred.second);
		}
	}
	
	void print_graph() {
		cout << fixed << setprecision(print_precision);
		cout << "\n\t\t\t Print Graph: \n";
		if (Graph.nodes.size() > 0) {
			cout << "\t\t\t\tnode\t\tedges\n";
			for (Dict<double,double>::Iterator _node = Graph.nodes.begin(); _node != Graph.nodes.end(); ++_node) {
				cout << "\t\t\t\t" << _node->first << " \t";
				dict_item<double, Dict<double,char>>* successors = Graph.succ.get(_node->second);
				if (successors != NULL) {
					Dict<double, char>* _succ = &successors->second;
					if (_succ->size() > 0) {
						for (auto& edges : *_succ)
							cout << edges.first << " " << edges.second << "; ";
					}
				}
				cout << "\n";
			}
		} else
			cout << "\t\t\t no graph to display\n\n";
		cout << "\n";
	}
	
	int size(){
		return Graph.nodes.size();
	}
	
	double start_node() {
		return Graph.nodes.min_entry();
	}
	
	double end_node() {
		return Graph.nodes.max_entry();
	}
};

#endif